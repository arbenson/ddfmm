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

#include "parUtils.h"

// TODO(arbenson): this is a bit of hack.
#define DVMAX 400

namespace par {
template <>
class Mpi_datatype<Index3> {
public:
    static MPI_Datatype value() {
        static bool         first = true;
        static MPI_Datatype datatype;
        if (first) {
            first = false;
            MPI_Type_contiguous(sizeof(Index3), MPI_BYTE, &datatype);
            MPI_Type_commit(&datatype);
        }
        return datatype;
    }
};
}


#ifdef LIMITED_MEMORY
bool CompareDownwardHighInfo(std::pair<double, Index3> a,
                             std::pair<double, Index3> b) {
    return a.first < b.first;
}

int Wave3d::LevelCommunication(std::map< double, std::vector<HFBoxAndDirectionKey> >& request_bnds,
                               double W) {
    std::vector<int> mask(HFBoxAndDirectionDat_Number, 0);
    mask[HFBoxAndDirectionDat_dirupeqnden] = 1;
    _bndvec.initialize_data();
    SAFE_FUNC_EVAL( _bndvec.getBegin(request_bnds[W], mask) );
    SAFE_FUNC_EVAL( _bndvec.getEnd(mask) );
    request_bnds[W].clear();
    std::ostringstream recv_msg;
    recv_msg << "kbytes received (W = " << W << ")";
    PrintCommData(GatherCommData(_bndvec.kbytes_received()), recv_msg.str());
    std::ostringstream sent_msg;
    sent_msg << "kbytes sent (W = " << W << ")";
    PrintCommData(GatherCommData(_bndvec.kbytes_sent()), sent_msg.str());
}
#endif

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

int Wave3d::LowFreqDownwardComm(std::set<BoxKey>& reqboxset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqDownwardComm");
#endif
    time_t t0 = time(0);
    std::vector<BoxKey> reqbox;
    reqbox.insert(reqbox.begin(), reqboxset.begin(), reqboxset.end());
    std::vector<int> mask(BoxDat_Number,0);
    mask[BoxDat_extden] = 1;
    mask[BoxDat_upeqnden] = 1;
    _boxvec.initialize_data();
    SAFE_FUNC_EVAL( _boxvec.getBegin(reqbox, mask) );
    SAFE_FUNC_EVAL( _boxvec.getEnd(mask) );
    time_t t1 = time(0);
    PrintParData(GatherParData(t0, t1), "Low frequency downward communication");
    PrintCommData(GatherCommData(_boxvec.kbytes_received()),
                  "kbytes received");
    PrintCommData(GatherCommData(_boxvec.kbytes_sent()),
                  "kbytes sent");
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

#ifdef LIMITED_MEMORY
int Wave3d::HighFreqPass(hdmap_t& hdmap) {
# ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqPass");
# endif
    time_t t0, t1;
    int mpirank = getMPIRank();
    
    if (mpirank == 0) {
        std::cout << "Beginning high frequency pass..." << std::endl;
    }
    t0 = time(0);

    // Find all directions on the first level (width = 1)
    std::vector<Index3> basedirs;
    double local_max_W = 1;
    for (hdmap_t::iterator mi = hdmap.begin(); mi != hdmap.end(); ++mi) {
        Index3 dir = mi->first;
        double W = dir2width(dir);
        if (W == 1) {
            basedirs.push_back(dir);
        }
        if (W > local_max_W) {
            local_max_W = W;
        }
    }

    double max_W = 1;
    MPI_Allreduce(&local_max_W, &max_W, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    std::set<HFBoxAndDirectionKey> reqbndset;
    for(int i = 0; i < basedirs.size(); ++i) {
        SAFE_FUNC_EVAL( EvalUpwardHighRecursive(1, basedirs[i], hdmap, reqbndset) );
    }
    t1 = time(0);
    PrintParData(GatherParData(t0, t1), "High frequency upward pass");

    t0 = time(0);
    std::vector<HFBoxAndDirectionKey> reqbnd;
    reqbnd.insert(reqbnd.begin(), reqbndset.begin(), reqbndset.end());

    // Requests by level
    std::map< double, std::vector<HFBoxAndDirectionKey> > request_bnds;
    while (!reqbnd.empty()) {
        HFBoxAndDirectionKey key = reqbnd.back();
        request_bnds[width(key.first)].push_back(key);
        reqbnd.pop_back();
    }

    std::list< std::vector< std::pair<double, Index3> > > all_info;
    for (int i = 0; i < basedirs.size(); ++i) {
        std::vector< std::pair<double, Index3> > curr_info;
        Index3 dir = basedirs[i];
        SAFE_FUNC_EVAL( GetDownwardHighInfo(1, dir, hdmap, curr_info) );
        std::sort(curr_info.begin(), curr_info.end(), &CompareDownwardHighInfo);
        all_info.push_back(curr_info);
    }

    // Level by level communication and computation
    for (double W = max_W; W >= 1; W /= 2) {
        LevelCommunication(request_bnds, W);
        time_t t2 = time(0);
        for (std::list< std::vector< std::pair<double, Index3> > >::iterator it = all_info.begin();
             it != all_info.end(); ++it) {
            while (!it->empty() && it->back().first == W) {
                Index3 dir = it->back().second;
                hdmap_t::iterator mi = hdmap.find(dir);
                CHECK_TRUE (mi != hdmap.end());
                EvalDownwardHigh(W, dir, mi->second);
                it->pop_back();
            }
        }
        time_t t3 =time(0);
        std::ostringstream msg;
        msg << "Computation time (W = " << W << ")";
        PrintParData(GatherParData(t2, t3), msg.str());
    }
    t1 = time(0);
    PrintParData(GatherParData(t0, t1), "High frequency downward pass");
    return 0;
}
#else
int Wave3d::HighFreqPass(hdmap_t& hdmap) {
# ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqPass");
# endif
    time_t t0, t1;
    int mpirank = getMPIRank();
    
    if(mpirank == 0) {
        std::cout << "Beginning high frequency pass..." << std::endl;
    }
    t0 = time(0);

    // Find all directions on the first level (width = 1)
    std::vector<Index3> basedirs;
    for (hdmap_t::iterator mi = hdmap.begin(); mi != hdmap.end(); ++mi) {
        Index3 dir = mi->first;
        if (Dir2Width(dir) == 1) {
            basedirs.push_back(dir);
        }
    }

    std::set<HFBoxAndDirectionKey> reqbndset;
    for (int i = 0; i < basedirs.size(); ++i) {
        SAFE_FUNC_EVAL( EvalUpwardHighRecursive(1, basedirs[i], hdmap, reqbndset) );
    }
    t1 = time(0);
    PrintParData(GatherParData(t0, t1), "High frequency upward pass");

    t0 = time(0);
    std::vector<int> mask(HFBoxAndDirectionDat_Number,0);
    mask[HFBoxAndDirectionDat_dirupeqnden] = 1;
    std::vector<HFBoxAndDirectionKey> reqbnd;
    reqbnd.insert(reqbnd.begin(), reqbndset.begin(), reqbndset.end());
    SAFE_FUNC_EVAL( _bndvec.getBegin(reqbnd, mask) );
    SAFE_FUNC_EVAL( _bndvec.getEnd(mask) );
    t1 = time(0);
    PrintParData(GatherParData(t0, t1), "High frequency communication");

    t0 = time(0);
    for (int i = 0; i < basedirs.size(); ++i) {
        Index3 dir = basedirs[i]; // LEXING: PRE HERE
        SAFE_FUNC_EVAL( EvalDownwardHighRecursive(1, dir, hdmap) );
    }
    t1 = time(0);

    PrintParData(GatherParData(t0, t1), "High frequency downward pass");
    return 0;
}
#endif

int Wave3d::GatherDensities(std::vector<int>& reqpts, ParVec<int,cpx,PtPrtn>& den) {
    int mpirank = getMPIRank();
    std::vector<int> all(1, 1);
    time_t t0 = time(0);
    SAFE_FUNC_EVAL( den.getBegin(reqpts, all) );
    SAFE_FUNC_EVAL( den.getEnd(all) );
    time_t t1 = time(0);
    if (mpirank == 0) {
        std::cout << "Density communication: " << difftime(t1, t0)
                  << " secs" << std::endl;
    }
    return 0;
}

int Wave3d::ConstructMaps(ldmap_t& ldmap, hdmap_t& hdmap) {
    int mpirank = getMPIRank();
    double eps = 1e-12;
    // construct maps, low frequency level by level, high frequency dir by dir
    // 
    // ldmap  maps box widths to a list of BoxKeys which correspond
    // to boxes in the low-frequency regime that are owned by this processor
    //
    // hdmap maps a direction to a pair of vectors:
    //     1. BoxKeys that have the direction in its outgoing direction list
    //     2. BoxKeys that have the direction in its incoming direction list
    //
    for (std::map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
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
                HFBoxAndDirectionDat dummy;
                // For each outgoing direction of this box, add to the first list
                for (std::set<Index3>::iterator si = curdat.outdirset().begin();
                    si != curdat.outdirset().end(); ++si) {
                    hdmap[*si].first.push_back(curkey);
                    // into bndvec
                    _bndvec.insert(HFBoxAndDirectionKey(curkey, *si), dummy);
                }
                // For each incoming direction of this box, add to the second list
                for (std::set<Index3>::iterator si = curdat.incdirset().begin();
                    si != curdat.incdirset().end(); ++si) {
                    hdmap[*si].second.push_back(curkey);
                    // into bndvec
                    _bndvec.insert(HFBoxAndDirectionKey(curkey, *si), dummy);
                }
            }
        }
    }
    return 0;
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
    std::vector<int> all(1, 1);
    ParVec<int, Point3, PtPrtn>& pos = (*_posptr);
    // Go through posptr to get nonlocal points
    std::vector<int> reqpts;
    for(std::map<int,Point3>::iterator mi = pos.lclmap().begin();
        mi != pos.lclmap().end(); ++mi) {
        reqpts.push_back( mi->first );
    }

    GatherDensities(reqpts, den);

    // Compute extden using ptidxvec
    for (std::map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (HasPoints(curdat) && OwnBox(curkey, mpirank) && IsTerminal(curdat)) {
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
    hdmap_t hdmap;
    ConstructMaps(ldmap, hdmap);

#if 0
    // Gather some statistics on the number of boxes and directions
    std::map<int, int> dircounts;
    for (std::map< Index3, box_lists_t >::iterator mi = hdmap.begin();
         mi != hdmap.end(); ++mi) {
        std::vector<BoxKey>& boxes = mi->second.first;
        for (int i = 0; i < boxes.size(); ++i) {
            dircounts[boxes[i].first] += BoxData(boxes[i]).DirInteractionListSize();
        } 
    }
    for (std::map<int, int>::iterator mi = dircounts.begin();
         mi != dircounts.end(); ++mi) {
        std::cerr << mpirank
                  << " directions: ("
                  << mi->first
                  << ", "
                  << mi->second
                  << ")"
                  << std::endl;
    }
#endif

    std::vector<Index3> box_locs;
    for (std::map<BoxKey, BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); ++mi) {
        int level = (mi->first).first;
        if (level == 5) {
            box_locs.push_back(mi->first.second);
        }
    }

    std::cerr << mpirank << ": first before sort: " << box_locs[0] << std::endl;
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    std::cerr << "Starting bitonicSort" << std::endl;
    par::bitonicSort(box_locs, MPI_COMM_WORLD);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    std::cerr << "Done with bitonicSort" << std::endl;
    std::cerr << mpirank << ": first after sort: " << box_locs[0] << std::endl;



    // Main work of the algorithm
    std::set<BoxKey> reqboxset;
    LowFreqUpwardPass(ldmap, reqboxset);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    HighFreqPass(hdmap);
    LowFreqDownwardComm(reqboxset);
    LowFreqDownwardPass(ldmap);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

    //set val from extval
    std::vector<int> wrtpts;
    for(std::map<int,Point3>::iterator mi = pos.lclmap().begin();
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
        if (HasPoints(curdat) && OwnBox(curkey, mpirank) && IsTerminal(curdat)) {
            CpxNumVec& extval = curdat.extval();
            std::vector<int>& curpis = curdat.ptidxvec();
            for (int k = 0; k < curpis.size(); ++k) {
                int poff = curpis[k];
                val.access(poff) = extval(k);
            }
        }
    }
    //call val->put
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
    //---------------
    int tdof = 1;
    for (int k = 0; k < srcvec.size(); ++k) {
        BoxKey srckey = srcvec[k];
        BoxDat& srcdat = _boxvec.access(srckey);
        CHECK_TRUE(HasPoints(srcdat));  // should have points

        Point3 srcctr = BoxCenter(srckey);
        //get array
        CpxNumVec upchkval(tdof*ucp.n());
        setvalue(upchkval,cpx(0,0));
        CpxNumVec& upeqnden = srcdat.upeqnden();
        //ue2dc
        if (IsTerminal(srcdat)) {
            DblNumMat upchkpos(ucp.m(), ucp.n());
            for (int k = 0; k < ucp.n(); ++k) {
                for (int d = 0; d < dim(); ++d) {
                    upchkpos(d,k) = ucp(d,k) + srcctr(d);
                }
            }
            //mul
            CpxNumMat mat;
            SAFE_FUNC_EVAL( _kernel.kernel(upchkpos, srcdat.extpos(), srcdat.extpos(), mat) );
            SAFE_FUNC_EVAL( zgemv(1.0, mat, srcdat.extden(), 1.0, upchkval) );
        } else {
            for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                int a = CHILD_IND1(ind);
                int b = CHILD_IND2(ind);
                int c = CHILD_IND3(ind);
                BoxKey key = ChildKey(srckey, Index3(a, b, c));
                std::pair<bool, BoxDat&> data = _boxvec.contains(key);
                if (data.first) {
                    BoxDat& chddat = data.second;
                    SAFE_FUNC_EVAL( zgemv(1.0, ue2uc(a, b, c), chddat.upeqnden(), 1.0, upchkval) );
                }
            }
        }

        //uc2ue
        CpxNumMat& v  = uc2ue(0);
        CpxNumMat& is = uc2ue(1); //LEXING: it is stored as a matrix
        CpxNumMat& up = uc2ue(2);
        CpxNumVec mid(up.m());
        setvalue(mid,cpx(0,0));
        SAFE_FUNC_EVAL( zgemv(1.0, up, upchkval, 0.0, mid) );
        for (int k = 0; k < mid.m(); ++k) {
            mid(k) = mid(k) * is(k,0);
        }
        upeqnden.resize(v.m());
        setvalue(upeqnden,cpx(0,0));
        SAFE_FUNC_EVAL( zgemv(1.0, v, mid, 0.0, upeqnden) );

        //-------------------------
        //EXTRA WORK, change role now
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
    int _P = P();
    for (int k = 0; k < trgvec.size(); ++k) {
        BoxKey trgkey = trgvec[k];
        BoxDat& trgdat = _boxvec.access(trgkey);
        CHECK_TRUE(HasPoints(trgdat));  // should have points

        Point3 trgctr = BoxCenter(trgkey);
        //array
        CpxNumVec& dnchkval = trgdat.dnchkval();
        if (dnchkval.m() == 0) {
            dnchkval.resize(dcp.n());
            setvalue(dnchkval,cpx(0,0));
        }
        if (trgdat.extval().m() == 0) {
            trgdat.extval().resize( trgdat.extpos().n() );
            setvalue(trgdat.extval(), cpx(0,0));
        }
        DblNumMat dnchkpos(dcp.m(), dcp.n());
        for (int k = 0; k < dcp.n(); ++k) {
            for (int d = 0; d < dim(); ++d) {
                dnchkpos(d,k) = dcp(d,k) + trgctr(d);
            }
        }
        // List computations
        SAFE_FUNC_EVAL( U_list_compute(trgdat) );
        SAFE_FUNC_EVAL( V_list_compute(trgdat, W, _P, trgctr, uep, dcp, dnchkval, ue2dc) );
        SAFE_FUNC_EVAL( W_list_compute(trgdat, W, uep) );
        SAFE_FUNC_EVAL( X_list_compute(trgdat, dcp, dnchkpos, dnchkval) );

        //-------------
        //dnchkval to dneqnden
        CpxNumMat& v  = dc2de(0);
        CpxNumMat& is = dc2de(1);
        CpxNumMat& up = dc2de(2);
        CpxNumVec mid(up.m());        setvalue(mid,cpx(0,0));
        SAFE_FUNC_EVAL( zgemv(1.0, up, dnchkval, 0.0, mid) );
        dnchkval.resize(0); //LEXING: SAVE SPACE
        for (int k = 0; k < mid.m(); ++k) {
            mid(k) = mid(k) * is(k,0);
        }
        CpxNumVec dneqnden(v.m());
        setvalue(dneqnden,cpx(0,0));
        SAFE_FUNC_EVAL( zgemv(1.0, v, mid, 0.0, dneqnden) );
        //-------------
        //to children or to exact points
        if (IsTerminal(trgdat)) {
            DblNumMat dneqnpos(dep.m(), dep.n());
            for (int k = 0; k < dep.n(); ++k) {
                for (int d = 0; d < dim(); ++d) {
                    dneqnpos(d,k) = dep(d,k) + trgctr(d);
                }
            }
            //mul
            CpxNumMat mat;
            SAFE_FUNC_EVAL( _kernel.kernel(trgdat.extpos(), dneqnpos, dneqnpos, mat) );
            SAFE_FUNC_EVAL( zgemv(1.0, mat, dneqnden, 1.0, trgdat.extval()) );
        } else {
            //put stuff to children
            for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                int a = CHILD_IND1(ind);
                int b = CHILD_IND2(ind);
                int c = CHILD_IND3(ind);
                BoxKey key = ChildKey(trgkey, Index3(a, b, c));
                std::pair<bool, BoxDat&> data = _boxvec.contains(key);
                if (!data.first) {
                  continue;
                }
                BoxDat& chddat = data.second;
                //mul
                if (chddat.dnchkval().m() == 0) {
                    chddat.dnchkval().resize(de2dc(a,b,c).m());
                    setvalue(chddat.dnchkval(), cpx(0,0));
                }
                SAFE_FUNC_EVAL( zgemv(1.0, de2dc(a, b, c), dneqnden, 1.0, chddat.dnchkval()) );
            }
        }
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalUpwardHighRecursive(double W, Index3 nowdir, hdmap_t& hdmap,
                                    std::set<HFBoxAndDirectionKey>& reqbndset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalUpwardHighRecursive");
#endif
    hdmap_t::iterator mi = hdmap.find(nowdir);
    if (mi != hdmap.end()) {
        SAFE_FUNC_EVAL( EvalUpwardHigh(W, nowdir, mi->second, reqbndset) );
        std::vector<Index3> dirvec = ChildDir(nowdir);
        for (int k = 0; k < dirvec.size(); ++k) {
            SAFE_FUNC_EVAL( EvalUpwardHighRecursive(2 * W, dirvec[k], hdmap, reqbndset) );
        }
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalDownwardHighRecursive(double W, Index3 nowdir, hdmap_t& hdmap) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalDownwardHighRecursive");
#endif
    hdmap_t::iterator mi = hdmap.find(nowdir);
    if (mi != hdmap.end()) {
        std::vector<Index3> dirvec = ChildDir(nowdir);
        for (int k = 0; k < dirvec.size(); ++k) {
            SAFE_FUNC_EVAL( EvalDownwardHighRecursive(2 * W, dirvec[k], hdmap) );
        }
        SAFE_FUNC_EVAL( EvalDownwardHigh(W, nowdir, mi->second) );
    }
    return 0;
}

//---------------------------------------------------------------------
# ifdef LIMITED_MEMORY
int Wave3d::GetDownwardHighInfo(double W, Index3 nowdir, hdmap_t& hdmap,
                                std::vector< std::pair<double, Index3> >& compute_info) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::DownwardHighCallStack");
#endif
    hdmap_t::iterator mi = hdmap.find(nowdir);
    if (mi != hdmap.end()) {
        std::vector<Index3> dirvec = chddir(nowdir);
        for (int k = 0; k < dirvec.size(); ++k) {
            SAFE_FUNC_EVAL( GetDownwardHighInfo(2 * W, dirvec[k], hdmap, compute_info) );
        }
        compute_info.push_back(std::pair<double, Index3>(W, nowdir));
    }
    return 0;
}
# endif

//---------------------------------------------------------------------
int Wave3d::EvalUpwardHigh(double W, Index3 dir, box_lists_t& hdvecs,
                           std::set<HFBoxAndDirectionKey>& reqbndset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalUpwardHigh");
#endif
    double eps = 1e-12;
    DblNumMat uep;
    DblNumMat ucp;
    NumVec<CpxNumMat> uc2ue;
    NumTns<CpxNumMat> ue2uc;
    SAFE_FUNC_EVAL( _mlibptr->UpwardHighFetch(W, dir, uep, ucp, uc2ue, ue2uc) );
    //---------------
    std::vector<BoxKey>& srcvec = hdvecs.first;
    for (int k = 0; k < srcvec.size(); ++k) {
        BoxKey srckey = srcvec[k];
        BoxDat& srcdat = _boxvec.access(srckey);
        CHECK_TRUE(HasPoints(srcdat));  // Should have points

        Point3 srcctr = BoxCenter(srckey);
        HFBoxAndDirectionKey bndkey(srckey, dir);
        HFBoxAndDirectionDat& bnddat = _bndvec.access( bndkey );
        CpxNumVec& upeqnden = bnddat.dirupeqnden();
        //eval
        CpxNumVec upchkval(ue2uc(0,0,0).m());
        setvalue(upchkval,cpx(0,0));
        // High-frequency M2M
        if (abs(W-1) < eps) {
            // The children boxes only have non-directional equivalent densities
            for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                int a = CHILD_IND1(ind);
                int b = CHILD_IND2(ind);
                int c = CHILD_IND3(ind);
                BoxKey key = ChildKey(srckey, Index3(a, b, c));
                // Do not compute unless _boxvec has the child key
                std::pair<bool, BoxDat&> data = _boxvec.contains(key);
                if (data.first) {
                    BoxDat& chddat = data.second;
                    CHECK_TRUE(HasPoints(chddat));
                    CpxNumVec& chdued = chddat.upeqnden();
                    SAFE_FUNC_EVAL( zgemv(1.0, ue2uc(a,b,c), chdued, 1.0, upchkval) );
                }
            }
        } else {
            // Pick the direction such that the child wedges in that direction
            // contain the parent wedge in direction dir
            Index3 pdir = ParentDir(dir);
            for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                int a = CHILD_IND1(ind);
                int b = CHILD_IND2(ind);
                int c = CHILD_IND3(ind);
                BoxKey chdkey = ChildKey(srckey, Index3(a, b, c));
                std::pair<bool, BoxDat&> data = _boxvec.contains(chdkey);
                if (data.first) {
                    BoxDat& chddat = data.second;
                    CHECK_TRUE(HasPoints(chddat));
                    HFBoxAndDirectionKey bndkey(chdkey, pdir);
                    HFBoxAndDirectionDat& bnddat = _bndvec.access(bndkey);
                    CpxNumVec& chdued = bnddat.dirupeqnden();
                    SAFE_FUNC_EVAL( zgemv(1.0, ue2uc(a,b,c), chdued, 1.0, upchkval) );
                }
            }
        }

        // Upward check to upward equivalency (uc2ue)
        CpxNumMat& E1 = uc2ue(0);
        CpxNumMat& E2 = uc2ue(1);
        CpxNumMat& E3 = uc2ue(2);
        cpx dat0[DVMAX], dat1[DVMAX];
        CpxNumVec tmp0(E3.m(), false, dat0);
        CHECK_TRUE(DVMAX >= E3.m());
        CpxNumVec tmp1(E2.m(), false, dat1);
        CHECK_TRUE(DVMAX >= E2.m());
        upeqnden.resize(E1.m());
        setvalue(upeqnden,cpx(0,0));
        SAFE_FUNC_EVAL( zgemv(1.0, E3, upchkval, 0.0, tmp0) );
        SAFE_FUNC_EVAL( zgemv(1.0, E2, tmp0, 0.0, tmp1) );
        SAFE_FUNC_EVAL( zgemv(1.0, E1, tmp1, 0.0, upeqnden) );
    }

    SAFE_FUNC_EVAL( GetInteractionListKeys(dir, hdvecs.second, reqbndset) );
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::GetInteractionListKeys(Index3 dir, std::vector<BoxKey>& target_boxes,
                                   std::set<HFBoxAndDirectionKey>& reqbndset) {
#ifndef RELEASE
  CallStackEntry entry("Wave3d::GetInteractionListKeys");
#endif
  for (int k = 0; k < target_boxes.size(); ++k) {
      BoxKey trgkey = target_boxes[k];
      BoxDat& trgdat = _boxvec.access(trgkey);
      CHECK_TRUE(HasPoints(trgdat));
      std::vector<BoxKey>& tmpvec = trgdat.fndeidxvec()[dir];
      for (int i = 0; i < tmpvec.size(); ++i) {
          BoxKey srckey = tmpvec[i];
          reqbndset.insert(HFBoxAndDirectionKey(srckey, dir));
      }
  }
  return 0;
}


//---------------------------------------------------------------------

//---------------------------------------------------------------------
int Wave3d::EvalDownwardHigh(double W, Index3 dir, box_lists_t& hdvecs) {
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
    std::vector<BoxKey>& trgvec = hdvecs.second;
    for (int k = 0; k < trgvec.size(); ++k) {
        BoxKey trgkey = trgvec[k];
        BoxDat& trgdat = _boxvec.access(trgkey);
        SAFE_FUNC_EVAL( HighFrequencyM2L(W, dir, trgkey, trgdat, dcp, uep) );
        SAFE_FUNC_EVAL( HighFrequencyL2L(W, dir, trgkey, dc2de, de2dc) );
     }

    // Now that we are done at this level, clear data to save on memory.
    std::vector<BoxKey>& srcvec = hdvecs.first;
    for (int k = 0; k < srcvec.size(); ++k) {
        BoxKey srckey = srcvec[k];
        BoxDat& srcdat = _boxvec.access(srckey);
        CHECK_TRUE(HasPoints(srcdat));  // should have points
        HFBoxAndDirectionKey bndkey(srckey, dir);
        HFBoxAndDirectionDat& bnddat = _bndvec.access( bndkey );
        bnddat.dirupeqnden().resize(0);
    }
    for (int k = 0; k < trgvec.size(); ++k) {
        BoxKey trgkey = trgvec[k];
        BoxDat& trgdat = _boxvec.access(trgkey);
        CHECK_TRUE(HasPoints(trgdat));  // should have points
        std::vector<BoxKey>& tmpvec = trgdat.fndeidxvec()[dir];
        for (int i = 0; i < tmpvec.size(); ++i) {
            BoxKey srckey = tmpvec[i];
            HFBoxAndDirectionKey bndkey(srckey, dir);
            HFBoxAndDirectionDat& bnddat = _bndvec.access(bndkey);
            bnddat.dirupeqnden().resize(0);
        }
    }
    return 0;
}
