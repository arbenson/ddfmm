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

int Wave3d::LowFreqUpwardPass(std::set<BoxKey>& reqboxset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqUpwardPass");
#endif
    int mpirank = getMPIRank();
    if (mpirank == 0) {
        std::cout << "Beginning low frequency upward pass..." << std::endl;
    }

    double t0 = MPI_Wtime();
    // For each box width in the low frequency regime that this processor
    // owns, evaluate upward.
    for (auto& kv : _ldmap) {
        SAFE_FUNC_EVAL( EvalUpwardLow(kv.first, kv.second, reqboxset) );
    }
    double t1 = MPI_Wtime();
    PrintParData(GatherParData(t0, t1), "Low frequency upward pass");
    return 0;
}

int Wave3d::LowFreqDownwardPass() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqDownwardPass");
#endif
    double t0 = MPI_Wtime();
    for (ldmap_t::reverse_iterator mi = _ldmap.rbegin();
        mi != _ldmap.rend(); ++mi) {
        SAFE_FUNC_EVAL( EvalDownwardLow(mi->first, mi->second) );
    }
    double t1 = MPI_Wtime();
    PrintParData(GatherParData(t0, t1), "Low frequency downward pass");
    return 0;
}

int Wave3d::HighFreqPass() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqPass");
#endif
    int mpirank = getMPIRank();

    level_hdkeys_map_t& level_hdmap_out = _level_prtns._level_hdmap_out;
    level_hdkeys_map_t& level_hdmap_inc = _level_prtns._level_hdmap_inc;
    
    if(mpirank == 0) {
        std::cout << "Beginning high frequency pass..." << std::endl;
    }

    // Upward pass (M2M)
    double t0, t1, t2, t3;
    t0 = MPI_Wtime();
    for (int level = level_hdmap_out.size() - 1; level >= 0; --level) {
        double W = _K / pow2(level);
        if (mpirank == 0) {
            std::cout << "---------------------------------" << std::endl;
            std::cout << "Box width for upwards pass: " << W << std::endl;
        }
        std::map<Index3, std::vector<BoxKey> >& level_out = level_hdmap_out[level];

	t2 = MPI_Wtime();
        for (auto& kv : level_out) {
            Index3 dir = kv.first;
            std::vector<BoxKey>& keys_out = kv.second;
            SAFE_FUNC_EVAL( EvalUpwardHigh(W, dir, keys_out) );
        }
	t3 = MPI_Wtime();
	PrintParData(GatherParData(t2, t3), "M2M level computation.");

        // Handle communication for this level.  We need to pass
        // directional upward equivalent densities from the children to the
        // parents on the next level. We assume that boxes on the unit level
        // are partitioned  by process, so we do not need to do any
        // communication for that level.
	t2 = MPI_Wtime();
        HighFreqM2MLevelComm(level);
	t3 = MPI_Wtime();
	PrintParData(GatherParData(t2, t3), "M2M level communication.");

        // TODO(arbenson): If we need memory, we can remove data received from other
	// processors at this level.
    }
    t1 = MPI_Wtime();
    PrintParData(GatherParData(t0, t1), "High frequency upward pass");

    // Communication for M2L
    t0 = MPI_Wtime();
    for (int level = 0;
         level < static_cast<int>(_level_prtns._hdkeys_inc.size()); ++level) {
        if (mpirank == 0) {
            std::cout << "----------------------" << std::endl;
            std::cout << "Level: " << level << std::endl;
	}
        std::set<BoxAndDirKey> request_keys;
        HighFreqInteractionListKeys(level, request_keys);
        HighFreqM2LComm(level, request_keys);
        SAFE_FUNC_EVAL(MPI_Barrier(MPI_COMM_WORLD));
    }
    t1 = MPI_Wtime();
    PrintParData(GatherParData(t0, t1), "High frequency M2L communication");

    // Downwards pass
    t0 = MPI_Wtime();
    for (int level = 0; level < static_cast<int>(level_hdmap_inc.size());
	 ++level) {
        double W = _K / pow2(level);
        SAFE_FUNC_EVAL(MPI_Barrier(MPI_COMM_WORLD));
        if (mpirank == 0) {
            std::cout << "---------------------------------" << std::endl;
            std::cout << "Box width for downwards pass: " << W << std::endl;
        }
        std::map<Index3, std::vector<BoxKey> >& level_inc = level_hdmap_inc[level];

        // Handle pre-communication for this level.  We send up the directional
        // downward check values for the children.  We assume that boxes on the
        // unit level are partitioned by process, so we do not need to do any
        // communication for that level.
	t2 = MPI_Wtime();
        HighFreqL2LLevelCommPre(level);
	t3 = MPI_Wtime();
	PrintParData(GatherParData(t2, t3), "L2L level communication (pre).");

        // Maintain set of keys that get updated by L2L
	t2 = MPI_Wtime();
        std::set<BoxAndDirKey> affected_keys;
        for (auto& kv : level_inc) {
            Index3 dir = kv.first;
            std::vector<BoxKey>& keys_inc = kv.second;
            SAFE_FUNC_EVAL( EvalDownwardHigh(W, dir, keys_inc, affected_keys) );
        }
	t3 = MPI_Wtime();
	PrintParData(GatherParData(t2, t3), "High frequency L2L compuation.");

        // Handle post-communication for this level.  We send back the
        // downward check values for the children. Again, we assume that the
        // boxes on the unit level are partitioned by process, so we do not
        // need to do any communication for that level.
	t2 = MPI_Wtime();
        HighFreqL2LLevelCommPost(level, affected_keys);
	t3 = MPI_Wtime();
	PrintParData(GatherParData(t2, t3), "L2L level communication (post).");

        // Remove old data from memory.
        CleanLevel(level);
    }
    t1 = MPI_Wtime();
    PrintParData(GatherParData(t0, t1), "High frequency downward pass");
    return 0;
}

int Wave3d::GatherLocalKeys() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::GatherLocalKeys");
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
    for (auto& kv : _boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        double W = BoxWidth(curkey);
        if (HasPoints(curdat) && OwnBox(curkey, mpirank)) {
            // Boxes of width less than one that are nonempty and are owned
            // by this processor get put in the low-frequency map.
          if (W > 1 - eps) {
                // High frequency regime
                BoxAndDirDat dummy;
                // For each outgoing direction of this box, add to the first list
                for (Index3 dir : curdat.outdirset()) {
                    // into bndvec
                    _bndvec.insert(BoxAndDirKey(curkey, dir), dummy);
                    int level = curkey._level;
                    level_hdkeys_out[level].push_back(BoxAndDirKey(curkey, dir));
                }
                
                // For each incoming direction of this box, add to the second list
                for (Index3 dir : curdat.incdirset()) {
                    BoxAndDirDat dat;
                    std::vector<BoxKey>& tmpvec = curdat.fndeidxvec()[dir];
                    std::vector<BoxAndDirKey>& interactionlist = dat.interactionlist();
                    interactionlist.resize(tmpvec.size());
                    for (int i = 0; i < static_cast<int>(tmpvec.size()); ++i) {
                        interactionlist[i] = BoxAndDirKey(tmpvec[i], dir);
                    }
                    // into bndvec
                    _bndvec.insert(BoxAndDirKey(curkey, dir), dat);
                    int level = curkey._level;
                    level_hdkeys_inc[level].push_back(BoxAndDirKey(curkey, dir));
                }

                // Save memory by clearing interaction lists stored in curdat.
                for (auto& kv : curdat.fndeidxvec()) {
                    std::vector<BoxKey>().swap(kv.second);
                }
            }
        }
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::eval(ParVec<int, cpx, PtPrtn>& den, ParVec<int, cpx, PtPrtn>& val) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::eval");
#endif
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    _self = this;
    int mpirank = getMPIRank();

    // Compute extden on leaf nodes using ptidxvec
    GatherDensities(den);

    // Main work of the algorithm
    double t0 = MPI_Wtime();
    std::set<BoxKey> reqboxset;
    LowFreqUpwardPass(reqboxset);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    HighFreqPass();
    LowFreqDownwardComm(reqboxset);
    LowFreqDownwardPass();
    double t1 = MPI_Wtime();
    PrintParData(GatherParData(t0, t1), "Actual wave evaluation computation time.");
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

    // For all values that I have but don't own, add to the list.
    std::vector<int> write_pts;
    for (auto& kv : _level_prtns._lf_boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        if (HasPoints(curdat) && 
	    _level_prtns.Owner(curkey) == mpirank &&
	    IsLeaf(curdat)) {
            CpxNumVec& extval = curdat.extval();
            std::vector<int>& curpis = curdat.ptidxvec();
            for (int k = 0; k < static_cast<int>(curpis.size()); ++k) {
                int poff = curpis[k];
                val.lclmap()[poff] = extval(k);
		if (val.prtn().owner(poff) != mpirank) {
                    write_pts.push_back(poff);
		}
            }
        }
    }
    std::vector<int> all(1, 1);
    val.putBegin(write_pts, all);
    val.putEnd(all);
    val.discard(write_pts);
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
    SAFE_FUNC_EVAL( _mlib.UpwardLowFetch(W, uep, ucp, uc2ue, ue2uc) );
    for (BoxKey& srckey : srcvec) {
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
    SAFE_FUNC_EVAL( _mlib.DownwardLowFetch(W, dep, dcp, dc2de, de2dc, ue2dc, uep) );
    //------------------
    for (BoxKey& trgkey : trgvec) {
        BoxDat& trgdat = _level_prtns._lf_boxvec.access(trgkey);
        CHECK_TRUE(HasPoints(trgdat));  // should have points
        CpxNumVec dneqnden;
        LowFreqM2L(trgkey, trgdat, dcp, ue2dc, dneqnden, uep, dc2de);
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
    SAFE_FUNC_EVAL( _mlib.UpwardHighFetch(W, dir, uep, ucp, uc2ue, ue2uc) );
    for (BoxKey& srckey : srcvec) {
        Point3 srcctr = BoxCenter(srckey);
        BoxAndDirKey bndkey(srckey, dir);
        HighFreqM2M(bndkey, uc2ue, ue2uc);
    }

    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalDownwardHigh(double W, Index3 dir, std::vector<BoxKey>& trgvec,
                             std::set<BoxAndDirKey>& affected_keys) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalDownwardHigh");
#endif
    DblNumMat dep;
    DblNumMat dcp;
    NumVec<CpxNumMat> dc2de;
    NumTns<CpxNumMat> de2dc;
    DblNumMat uep;
    SAFE_FUNC_EVAL( _mlib.DownwardHighFetch(W, dir, dep, dcp, dc2de, de2dc, uep) );
    for (BoxKey& trgkey : trgvec) {
        BoxAndDirKey bndkey(trgkey, dir);
        SAFE_FUNC_EVAL( HighFreqM2L(bndkey, dcp, uep) );
        SAFE_FUNC_EVAL( HighFreqL2L(bndkey, dc2de, de2dc, affected_keys) );
    }
    return 0;
}
