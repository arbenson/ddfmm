/* Distributed Directional Fast Multipole Method
   Copyright (C) 2014 Austin Benson, Lexing Ying, and Jack Poulson

 This file is part of DDFMM.

    DDFMM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DDFMM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DDFMM.  If not, see <http://www.gnu.org/licenses/>. */
#include "wave3d.hpp"
#include "vecmatop.hpp"
#include "DataCollection.hpp"

#include <algorithm>
#include <list>
#include <map>
#include <sstream>
#include <vector>

// TODO(arbenson): this is a bit of hack.
#define DVMAX 400

int Wave3d::LowFreqUpwardPass(ldmap_t& ldmap, std::set<BoxKey>& reqboxset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqUpwardPass");
#endif
    int mpirank = getMPIRank();
    if (mpirank == 0) {
        std::cout << "Beginning low frequency upward pass..." << std::endl;
    }

    time_t t0 = time(0);
    // For each box width in the low frequency regime that this processor
    // owns, evaluate upward.
    for (ldmap_t::iterator mi = ldmap.begin(); mi != ldmap.end(); ++mi) {
        SAFE_FUNC_EVAL( EvalUpwardLow(mi->first, mi->second, reqboxset) );
    }
    time_t t1 = time(0);
    PrintParData(GatherParData(t0, t1), "Low frequency upward pass");
    return 0;
}

int Wave3d::LowFreqDownwardPass(ldmap_t& ldmap) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqDownwardPass");
#endif
    time_t t0 = time(0);
    for (ldmap_t::reverse_iterator mi = ldmap.rbegin();
        mi != ldmap.rend(); ++mi) {
        SAFE_FUNC_EVAL( EvalDownwardLow(mi->first, mi->second) );
    }
    time_t t1 = time(0);
    PrintParData(GatherParData(t0, t1), "Low frequency downward pass");
    return 0;
}

int Wave3d::HighFreqPass(level_hdkeys_map_t& level_hdmap_out,
                         level_hdkeys_map_t& level_hdmap_inc) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqPass");
#endif
    time_t t0, t1;
    int mpirank = getMPIRank();
    
    if(mpirank == 0) {
        std::cout << "Beginning high frequency pass..." << std::endl;
    }

    // Upward pass (M2M)
    t0 = time(0);

    for (int level = level_hdmap_out.size() - 1; level >= 0; --level) {
        double W = _K / pow2(level);
        std::map<Index3, std::vector<BoxKey> >& level_out = level_hdmap_out[level];

        // Handle communication for this level.  We need request the directional
        // upward equivalent densities from the children needed on this level.
        // We assume that boxes on the unit level are partitioned by process,
        // so we do not need to do any communication for that level.
        if (level < UnitLevel()) {
            HighFreqM2MLevelComm(level);
        }

        for (std::map<Index3, std::vector<BoxKey> >::iterator mi = level_out.begin();
            mi != level_out.end(); ++mi) {
            Index3 dir = mi->first;
            std::vector<BoxKey>& keys_out = mi->second;
            SAFE_FUNC_EVAL( EvalUpwardHigh(W, dir, keys_out) );
        }

        // TODO(arbenson): remove data from parvec to save memory
    }
    t1 = time(0);
    PrintParData(GatherParData(t0, t1), "High frequency upward pass");

    // Communication for M2L
    t0 = time(0);
    for (int level = 0; level < _level_prtns._hdkeys_inc.size(); ++level) {
        std::set<BoxAndDirKey> request_keys;
        HighFreqInteractionListKeys(level, request_keys);
	HighFreqM2LComm(level, request_keys);
        SAFE_FUNC_EVAL(MPI_Barrier(MPI_COMM_WORLD));
    }
    t1 = time(0);
    PrintParData(GatherParData(t0, t1), "High frequency communication");

    // Downwards pass
    t0 = time(0);
    for (int level = 0; level < level_hdmap_inc.size(); ++level) {
        double W = _K / pow2(level);
        std::map<Index3, std::vector<BoxKey> > level_inc = level_hdmap_inc[level];
        std::map<Index3, std::vector<BoxKey> > level_out = level_hdmap_out[level];

        // Handle pre-communication for this level.  We request the directional
        // downward check values for the children.
        // We assume that boxes on the unit level are partitioned by process,
        // so we do not need to do any communication for that level.
        if (level < UnitLevel()) {
            HighFreqL2LLevelCommPre(level);
        }

        for (std::map<Index3, std::vector<BoxKey> >::iterator mi = level_inc.begin();
            mi != level_inc.end(); ++mi) {
            Index3 dir = mi->first;
            std::vector<BoxKey>& keys_inc = mi->second;
            std::vector<BoxKey>& keys_out = level_out[dir];
            SAFE_FUNC_EVAL( EvalDownwardHigh(W, dir, keys_inc, keys_out) );
        }

        // Handle post-communication for this level.  We send back the
        // downward check values for the children.
        // We assume that boxes on the unit level are partitioned by process,
        // so we do not need to do any communication for that level.
        if (level < UnitLevel()) {
            HighFreqL2LLevelCommPost(level);
        }

        // TODO(arbenson): remove data from parvec to save memory
    }
    t1 = time(0);

    PrintParData(GatherParData(t0, t1), "High frequency downward pass");
    return 0;
}

int Wave3d::ConstructMaps(ldmap_t& ldmap,
                          level_hdkeys_map_t& level_hdmap_out,
                          level_hdkeys_map_t& level_hdmap_inc) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::ConstructMaps");
#endif
    int mpirank = getMPIRank();
    double eps = 1e-12;
    level_hdkeys_t& level_hdkeys_out = _level_prtns._hdkeys_out;
    level_hdkeys_t& level_hdkeys_inc = _level_prtns._hdkeys_inc;

    // construct maps, low frequency level by level, high frequency dir by dir
    // 
    // ldmap  maps box widths to a list of BoxKeys which correspond
    // to boxes in the low-frequency regime that are owned by this processor
    //
    // hdmap maps a direction to a pair of vectors:
    //     1. BoxKeys that have the direction in its outgoing direction list
    //     2. BoxKeys that have the direction in its incoming direction list
    //
    for (std::map<BoxKey, BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        double W = BoxWidth(curkey);
        if (HasPoints(curdat) && OwnBox(curkey, mpirank)) {
            // Boxes of width less than one that are nonempty and are owned
            // by this processor get put in the low-frequency map.
            if (W < 1 - eps) {
                ldmap[W].push_back(curkey);
	    } else {
                // High frequency regime
                BoxAndDirDat dummy;
                // For each outgoing direction of this box, add to the first list
                for (std::set<Index3>::iterator si = curdat.outdirset().begin();
                    si != curdat.outdirset().end(); ++si) {
                    // into bndvec
                    _bndvec.insert(BoxAndDirKey(curkey, *si), dummy);
                    int level = curkey.first;
                    level_hdkeys_out[level].push_back(BoxAndDirKey(curkey, *si));
                    level_hdmap_out[level][*si].push_back(curkey);
                }
                
                // For each incoming direction of this box, add to the second list
                for (std::set<Index3>::iterator si = curdat.incdirset().begin();
                    si != curdat.incdirset().end(); ++si) {
                    BoxAndDirDat dat;
                    Index3 dir = *si;
                    std::vector<BoxKey>& tmpvec = curdat.fndeidxvec()[dir];
                    std::vector<BoxAndDirKey>& interactionlist = dat.interactionlist();
                    interactionlist.resize(tmpvec.size());
                    for (int i = 0; i < tmpvec.size(); ++i) {
                        interactionlist[i] = BoxAndDirKey(tmpvec[i], dir);
                    }
                    // into bndvec
                    _bndvec.insert(BoxAndDirKey(curkey, dir), dat);
                    int level = curkey.first;
                    level_hdkeys_inc[level].push_back(BoxAndDirKey(curkey, *si));
                    level_hdmap_inc[level][dir].push_back(curkey);
                }

                // Save memory by clearing interaction lists stored in curdat.
                for (std::map< Index3, std::vector<BoxKey> >::iterator mi = curdat.fndeidxvec().begin();
                     mi != curdat.fndeidxvec().end(); ++mi) {
                    std::vector<BoxKey>& tmpvec = mi->second;
                    std::vector<BoxKey>().swap(tmpvec);
                }
            }
        }
    }
    return 0;
}

int Wave3d::ConstructMaps2(ldmap_t& ldmap) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::ConstructMaps2");
#endif
    ldmap.clear();
    int mpirank = getMPIRank();
    double eps = 1e-12;
    for (std::map<BoxKey, BoxDat>::iterator mi = _level_prtns._lf_boxvec.lclmap().begin();
	 mi != _level_prtns._lf_boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        double W = BoxWidth(curkey);
        if (HasPoints(curdat) && _level_prtns._lf_boxvec.prtn().owner(curkey) == mpirank) {
	    CHECK_TRUE(W < 1 - eps);
            ldmap[W].push_back(curkey);
	}
    }
}


//---------------------------------------------------------------------
int Wave3d::eval(ParVec<int,cpx,PtPrtn>& den, ParVec<int,cpx,PtPrtn>& val) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::eval");
#endif
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    _self = this;
    time_t t0, t1, t2, t3;
    int mpirank = getMPIRank();

#if 0
    GatherDensities(reqpts, den);

    // Compute extden on leaf nodes using ptidxvec
    for (std::map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (HasPoints(curdat) && OwnBox(curkey, mpirank) && IsLeaf(curdat)) {
            std::vector<int>& curpis = curdat.ptidxvec();
            CpxNumVec& extden = curdat.extden();
            extden.resize(curpis.size());
            for (int k = 0; k < curpis.size(); ++k) {
                int poff = curpis[k];
                extden(k) = den.access(poff);
            }
        }
    }
    SAFE_FUNC_EVAL( den.discard(reqpts) );
#endif

    // Delete of empty boxes
    std::list<BoxKey> to_delete;
    for (std::map<BoxKey, BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (!HasPoints(curdat)) {
            to_delete.push_back(curkey);
        }
    }
    std::cerr << "Deleting " << to_delete.size()
              << " out of " << _boxvec.lclmap().size()
              << " total boxes." << std::endl;
    for (std::list<BoxKey>::iterator mi = to_delete.begin();
         mi != to_delete.end(); ++mi) {
        BoxKey curkey = *mi;
        _boxvec._lclmap.erase(curkey);
    }

    // Setup of low and high frequency maps
    ldmap_t ldmap;
    int max_level = 10;
    _level_prtns.init(max_level);
    level_hdkeys_map_t level_hdmap_out(max_level);
    level_hdkeys_map_t level_hdmap_inc(max_level);
    ConstructMaps(ldmap, level_hdmap_out, level_hdmap_inc);

    PrtnDirections(_level_prtns._hdkeys_out,
                   _level_prtns._hf_vecs_out);
    PrtnDirections(_level_prtns._hdkeys_inc,
                   _level_prtns._hf_vecs_inc);
    PrtnUnitLevel();
    _level_prtns.FormMaps();

    // Gather box data at the unit level for the partitioning of the trees.
    std::vector<int> mask1(BoxDat_Number, 0);
    mask1[BoxDat_tag] = 1;
    mask1[BoxDat_ptidxvec] = 1;
    SAFE_FUNC_EVAL(_boxvec.getBegin(&Wave3d::TransferUnitLevelData_wrapper, mask1));
    SAFE_FUNC_EVAL(_boxvec.getEnd(mask1));

    std::vector<int> mask2(BoxAndDirDat_Number, 0);
    mask2[BoxAndDirDat_interactionlist] = 1;
    SAFE_FUNC_EVAL(_bndvec.getBegin(&Wave3d::TransferBoxAndDirData_wrapper, mask2));
    SAFE_FUNC_EVAL(_bndvec.getEnd(mask2));
    TransferDataToLevels();

    // Now we have the unit level information, so we can setup our part of the tree.
    SetupLowFreqOctree();

    // Compute extden on leaf nodes using ptidxvec
    GatherDensities2(den);
    for (std::map<BoxKey,BoxDat>::iterator mi = _level_prtns._lf_boxvec.lclmap().begin();
        mi != _level_prtns._lf_boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (HasPoints(curdat) &&
	    _level_prtns._lf_boxvec.prtn().owner(curkey) == mpirank &&
	    IsLeaf(curdat)) {
            std::vector<int>& curpis = curdat.ptidxvec();
            CpxNumVec& extden = curdat.extden();
            extden.resize(curpis.size());
            for (int k = 0; k < curpis.size(); ++k) {
                int poff = curpis[k];
                extden(k) = den.access(poff);
            }
        }
    }
#if 0
    // TODO(arbenson): discard these
    SAFE_FUNC_EVAL( den.discard(reqpts) );
#endif

    ConstructMaps2(ldmap);

    // Main work of the algorithm
    std::set<BoxKey> reqboxset;
    LowFreqUpwardPass(ldmap, reqboxset);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    HighFreqPass(_level_prtns._level_hdmap_out, _level_prtns._level_hdmap_inc);
    LowFreqDownwardComm(reqboxset);
    LowFreqDownwardPass(ldmap);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

    ParVec<int, Point3, PtPrtn>& pos = (*_posptr);
    std::vector<int> all(1, 1);
    //set val from extval
    std::vector<int> wrtpts;
    for(std::map<int, Point3>::iterator mi = pos.lclmap().begin();
        mi != pos.lclmap().end(); ++mi) {
        if (pos.prtn().owner(mi->first) != mpirank) {
            wrtpts.push_back(mi->first);
        }
    }
    val.expand(wrtpts);
    for (std::map<BoxKey, BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (HasPoints(curdat) && OwnBox(curkey, mpirank) && IsLeaf(curdat)) {
            CpxNumVec& extval = curdat.extval();
            std::vector<int>& curpis = curdat.ptidxvec();
            for (int k = 0; k < curpis.size(); ++k) {
                int poff = curpis[k];
                val.access(poff) = extval(k);
            }
        }
    }
    val.putBegin(wrtpts, all);  val.putEnd(all);
    val.discard(wrtpts);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalUpwardLow(double W, std::vector<BoxKey>& srcvec,
                          std::set<BoxKey>& reqboxset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalUpwardLow");
#endif
    DblNumMat uep;
    DblNumMat ucp;
    NumVec<CpxNumMat> uc2ue;
    NumTns<CpxNumMat> ue2uc;
    SAFE_FUNC_EVAL( _mlibptr->UpwardLowFetch(W, uep, ucp, uc2ue, ue2uc) );
    for (int k = 0; k < srcvec.size(); ++k) {
        BoxKey srckey = srcvec[k];
        BoxDat& srcdat = _level_prtns._lf_boxvec.access(srckey);
        LowFreqM2M(srckey, srcdat, uep, ucp, uc2ue, ue2uc);

        // Add boxes in U, V, W, and X lists of trgdat to reqboxset.           
        BoxDat& trgdat = srcdat;
        reqboxset.insert(trgdat.undeidxvec().begin(), trgdat.undeidxvec().end());
        reqboxset.insert(trgdat.vndeidxvec().begin(), trgdat.vndeidxvec().end());
        reqboxset.insert(trgdat.wndeidxvec().begin(), trgdat.wndeidxvec().end());
        reqboxset.insert(trgdat.xndeidxvec().begin(), trgdat.xndeidxvec().end());
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalDownwardLow(double W, std::vector<BoxKey>& trgvec) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalDownwardLow");
#endif
    DblNumMat dep;
    DblNumMat dcp;
    NumVec<CpxNumMat> dc2de;
    NumTns<CpxNumMat> de2dc;
    NumTns<CpxNumTns> ue2dc;
    DblNumMat uep;
    SAFE_FUNC_EVAL( _mlibptr->DownwardLowFetch(W, dep, dcp, dc2de, de2dc, ue2dc, uep) );
    //------------------
    for (int k = 0; k < trgvec.size(); ++k) {
        BoxKey trgkey = trgvec[k];
        BoxDat& trgdat = _level_prtns._lf_boxvec.access(trgkey);
        CHECK_TRUE(HasPoints(trgdat));  // should have points
        CpxNumVec dneqnden;
        LowFreqM2L(W, trgkey, trgdat, dcp, ue2dc, dneqnden, uep, dc2de);
        LowFreqL2L(trgkey, trgdat, dep, de2dc, dneqnden);
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalUpwardHigh(double W, Index3 dir, std::vector<BoxKey>& srcvec) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalUpwardHigh");
#endif
    DblNumMat uep;
    DblNumMat ucp;
    NumVec<CpxNumMat> uc2ue;
    NumTns<CpxNumMat> ue2uc;
    SAFE_FUNC_EVAL( _mlibptr->UpwardHighFetch(W, dir, uep, ucp, uc2ue, ue2uc) );
    for (int k = 0; k < srcvec.size(); ++k) {
        BoxKey srckey = srcvec[k];
        BoxDat& srcdat = _boxvec.access(srckey);
        CHECK_TRUE(HasPoints(srcdat));  // Should have points

        Point3 srcctr = BoxCenter(srckey);
        BoxAndDirKey bndkey(srckey, dir);
        HighFreqM2M(W, bndkey, uc2ue, ue2uc);
    }

    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalDownwardHigh(double W, Index3 dir, std::vector<BoxKey>& trgvec,
                             std::vector<BoxKey>& srcvec) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalDownwardHigh");
#endif
    int mpirank = getMPIRank();

    DblNumMat dep;
    DblNumMat dcp;
    NumVec<CpxNumMat> dc2de;
    NumTns<CpxNumMat> de2dc;
    DblNumMat uep;
    SAFE_FUNC_EVAL( _mlibptr->DownwardHighFetch(W, dir, dep, dcp, dc2de, de2dc, uep) );
    //LEXING: IMPORTANT
    for (int k = 0; k < trgvec.size(); ++k) {
        BoxKey trgkey = trgvec[k];
        BoxDat& trgdat = _boxvec.access(trgkey);
        SAFE_FUNC_EVAL( HighFreqM2L(W, dir, trgkey, trgdat, dcp, uep) );
        SAFE_FUNC_EVAL( HighFreqL2L(W, dir, trgkey, dc2de, de2dc) );
     }

    // Now that we are done at this level, clear data to save on memory.
    for (int k = 0; k < srcvec.size(); ++k) {
        BoxKey srckey = srcvec[k];
        BoxDat& srcdat = _boxvec.access(srckey);
        CHECK_TRUE(HasPoints(srcdat));  // should have points
        BoxAndDirKey bndkey(srckey, dir);
        BoxAndDirDat& bnddat = _bndvec.access( bndkey );
        bnddat.dirupeqnden().resize(0);
    }
    for (int k = 0; k < trgvec.size(); ++k) {
        BoxKey trgkey = trgvec[k];
        BoxDat& trgdat = _boxvec.access(trgkey);
        CHECK_TRUE(HasPoints(trgdat));  // should have points
        std::vector<BoxKey>& tmpvec = trgdat.fndeidxvec()[dir];
        for (int i = 0; i < tmpvec.size(); ++i) {
            BoxKey srckey = tmpvec[i];
            BoxAndDirKey bndkey(srckey, dir);
            BoxAndDirDat& bnddat = _bndvec.access(bndkey);
            bnddat.dirupeqnden().resize(0);
        }
    }
    return 0;
}
