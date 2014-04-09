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

//---------------------------------------------------------------------
int Wave3d::setup(std::map<std::string, std::string>& opts) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup");
#endif
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    _self = this;
    int mpirank = getMPIRank();
    // Read optional data
    std::map<std::string, std::string>::iterator mi;
    mi = opts.find("-" + prefix() + "ACCU");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _ACCU;
    }
    mi = opts.find("-" + prefix() + "NPQ");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _NPQ;
    }
    mi = opts.find("-" + prefix() + "K");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _K;
    }
    mi = opts.find("-" + prefix() + "ctr");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        double x, y, z;
        ss >> x >> y >> z;
        _ctr = Point3(x,y,z);
    }
    mi = opts.find("-" + prefix() + "ptsmax");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _ptsmax;
    }
    mi = opts.find("-" + prefix() + "maxlevel");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _maxlevel;
    }
    //
    if (mpirank == 0) {
        std::cout << _K <<      " | "
                  << _ACCU <<   " | "
                  << _NPQ <<    " | "
                  << _ctr <<    " | "
                  << _ptsmax << " | "
                  << _maxlevel
                  << std::endl;
    }
    //
    //create the parvecs
    BoxPrtn bp;
    bp.ownerinfo() = _geomprtn;
    _boxvec.prtn() = bp;
    BoxAndDirPrtn tp;
    tp.ownerinfo() = _geomprtn;
    _bndvec.prtn() = tp;
    //generate octree
    SAFE_FUNC_EVAL(SetupTree());
    //plans
    int _P = P();
    _denfft.resize(2*_P, 2*_P, 2*_P);
    _fplan = fftw_plan_dft_3d(2 * _P, 2 * _P, 2 * _P, (fftw_complex*) (_denfft.data()),
                              (fftw_complex*)(_denfft.data()), FFTW_FORWARD,
                              FFTW_MEASURE);
    CHECK_TRUE(_fplan != NULL);
    setvalue(_denfft,cpx(0,0));
    //
    _valfft.resize(2 * _P, 2 * _P, 2 * _P);
    _bplan = fftw_plan_dft_3d(2 * _P, 2 * _P, 2 * _P,
                              (fftw_complex*) (_valfft.data()),
                              (fftw_complex*) (_valfft.data()), FFTW_BACKWARD,
                                               FFTW_ESTIMATE); 
    CHECK_TRUE(_bplan != NULL);
    setvalue(_valfft,cpx(0,0));
    return 0;
}

int Wave3d::RecursiveBoxInsert(std::queue< std::pair<BoxKey, BoxDat> >& tmpq) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::RecursiveBoxInsert");
#endif
    ParVec<int, Point3, PtPrtn>& pos = *_posptr;
    while (!tmpq.empty()) {
        std::pair<BoxKey, BoxDat> curent = tmpq.front();
        tmpq.pop();
        BoxKey& curkey = curent.first;
        BoxDat& curdat = curent.second;
        //LEXING: VERY IMPORTANT
        if (curdat.ptidxvec().size() > 0) {
            curdat.tag() |=  WAVE3D_PTS;
        }
        // We take an action if we are in the high frequency regime with points OR
        // we are in the low frequency regime with sufficient number of points and
	// not at the max depth.
	// If no action is taken, we are at a leaf node.
        bool action = (curkey.first <= UnitLevel() && curdat.ptidxvec().size() > 0) ||
            (curdat.ptidxvec().size() > ptsmax() && curkey.first < maxlevel() - 1);
        if (action) {
            // Subdivide to get new children
            NumTns<BoxDat> chdboxtns(2, 2, 2);
            Point3 curctr = BoxCenter(curkey); //LEXING: VERY IMPORTANT
            for (int g = 0; g < curdat.ptidxvec().size(); ++g) {
                int tmpidx = curdat.ptidxvec()[g];
                Point3 tmp = pos.access(tmpidx); //get position value
                Index3 idx;
                for (int d = 0; d < 3; ++d) {
                    idx(d) = (tmp(d) >= curctr(d));
                }
                // put points to children
                chdboxtns(idx(0), idx(1), idx(2)).ptidxvec().push_back(tmpidx);
            }
            // Put non-empty ones into queue
            for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                int a = CHILD_IND1(ind);
                int b = CHILD_IND2(ind);
                int c = CHILD_IND3(ind);
                BoxKey key = ChildKey(curkey, Index3(a,b,c));
		tmpq.push( std::pair<BoxKey, BoxDat>(key, chdboxtns(a, b, c)) );
            }
            // Destory ptidxvec to save memory.
	    std::vector<int>().swap(curdat.ptidxvec());
        } else {
            // Copy data into _extpos
            curdat.extpos().resize(3, curdat.ptidxvec().size());
            for (int g = 0; g < curdat.ptidxvec().size(); ++g) {
                int tmpidx = curdat.ptidxvec()[g];
                Point3 tmp = pos.access(tmpidx);
                for (int d = 0; d < 3; ++d) {
                    curdat.extpos()(d, g) = tmp(d);
                }
            }
            //LEXING: VERY IMPORTANT
            curdat.tag() |= WAVE3D_LEAF;
        }
        // Add my self into _tree
        _boxvec.insert(curkey, curdat); //LEXING: CHECK
    }
}

//---------------------------------------------------------------------
int Wave3d::SetupTree() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupTree");
#endif
    int mpirank = getMPIRank();
    double eps = 1e-12;
    double K = this->K();
    // pos contains all data read by this processor
    ParVec<int, Point3, PtPrtn>& pos = *_posptr;

    // Get all of the geometry information needed for this processor
    std::vector<int> all(1, 1);
    SAFE_FUNC_EVAL( pos.getBegin(&(Wave3d::setup_Q1_wrapper), all) );
    SAFE_FUNC_EVAL( pos.getEnd(all) );

    int numC = _geomprtn.m();
    int lvlC = CellLevel();
    // Generate cell level boxes, put them into queue
    Point3 bctr = ctr();  // overall center of domain
    NumTns<BoxDat> cellboxtns(numC, numC, numC);
    // Fill boxes with points.
    for (std::map<int,Point3>::iterator mi = pos.lclmap().begin();
         mi != pos.lclmap().end(); ++mi) {
        int key = mi->first;
        Point3 pos = mi->second;
        Index3 idx;
        for(int d = 0; d < 3; d++) {
            idx(d) = (int) floor(numC * ((pos(d) - bctr(d) + K / 2) / K));
            CHECK_TRUE(idx(d) >= 0 && idx(d) < numC);
        }
        // Put the points in
        cellboxtns(idx(0), idx(1), idx(2)).ptidxvec().push_back( key );
    }
    // Put all boxes owned by this process (whether or not is empty) in a queue.
    // TODO(arbenson): this should be more efficient.
    std::queue< std::pair<BoxKey,BoxDat> > tmpq;
    for (int a = 0; a < numC; ++a) {
        for (int b = 0; b < numC; ++b) {
            for (int c = 0; c < numC; ++c) {
                if (_geomprtn(a, b, c) == mpirank) {
                    BoxKey key(lvlC, Index3(a, b, c));
                    tmpq.push( std::pair<BoxKey,BoxDat>(key, cellboxtns(a, b, c)) );
                }
            }
        }
    }
    cellboxtns.resize(0, 0, 0);

    // Construct the tree.
    RecursiveBoxInsert(tmpq);
 
    // call get setup_Q2
    std::vector<int> mask1(BoxDat_Number,0);
    mask1[BoxDat_tag] = 1;
    SAFE_FUNC_EVAL( _boxvec.getBegin( &(Wave3d::setup_Q2_wrapper), mask1 ) );
    SAFE_FUNC_EVAL( _boxvec.getEnd( mask1 ) );
    // Compute lists, low list and high list
    for (std::map<BoxKey, BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); mi++) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        // For all of my boxes with points, setup the call list.
        if (OwnBox(curkey, mpirank) && HasPoints(curdat)) {
            if (BoxWidth(curkey) < 1 - eps) { // strictly < 1
                // Low frequency regime
                SAFE_FUNC_EVAL(SetupTreeLowFreqLists(curkey, curdat));
            } else {
                // High frequency regime
                SAFE_FUNC_EVAL(SetupTreeHighFreqLists(curkey, curdat));
            }
        }
    }

    // Delete endeidxvec since it was only used to build the interaction lists
    // in the high-frequency regime.
    for (std::map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (OwnBox(curkey, mpirank)) {
            // Save memory
            std::vector<BoxKey>().swap(curdat.endeidxvec());
        }
    }

    // 3. get extpos
    std::set<BoxKey> reqboxset;
    for (std::map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (HasPoints(curdat) && OwnBox(curkey, mpirank)) {
            reqboxset.insert(curdat.undeidxvec().begin(), curdat.undeidxvec().end());
            reqboxset.insert(curdat.vndeidxvec().begin(), curdat.vndeidxvec().end());
            reqboxset.insert(curdat.wndeidxvec().begin(), curdat.wndeidxvec().end());
            reqboxset.insert(curdat.xndeidxvec().begin(), curdat.xndeidxvec().end());
        }
    }
    std::vector<BoxKey> reqbox;
    reqbox.insert(reqbox.begin(), reqboxset.begin(), reqboxset.end());
    std::vector<int> mask2(BoxDat_Number,0);
    mask2[BoxDat_extpos] = 1;
    SAFE_FUNC_EVAL( _boxvec.getBegin(reqbox, mask2) );
    SAFE_FUNC_EVAL( _boxvec.getEnd(mask2) );
    for (std::map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (HasPoints(curdat) && OwnBox(curkey, mpirank)) {
           for (std::vector<BoxKey>::iterator vi = curdat.vndeidxvec().begin();
                vi != curdat.vndeidxvec().end(); ++vi) {
                BoxKey neikey = (*vi);
                BoxDat& neidat = _boxvec.access(neikey);
                neidat.fftnum() ++;
            }
        }
    }

    // 4. dirupeqndenvec, dirdnchkvalvec
    for (std::map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); mi++) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        double W = BoxWidth(curkey);
        if (OwnBox(curkey, mpirank) && W > 1 - eps && HasPoints(curdat)) {
            if (!IsCellLevelBox(curkey)) {
                BoxKey parkey = ParentKey(curkey);
                BoxDat& pardat = BoxData(parkey);
                for (std::set<Index3>::iterator si = pardat.outdirset().begin();
                    si != pardat.outdirset().end(); si++) {
                    Index3 nowdir = ParentDir(*si);
                    curdat.outdirset().insert(nowdir);
                }
                for (std::set<Index3>::iterator si = pardat.incdirset().begin();
                    si != pardat.incdirset().end(); si++) {
                    Index3 nowdir = ParentDir(*si);
                    curdat.incdirset().insert(nowdir);
                }
            }
            Point3 curctr = BoxCenter(curkey);
            for (std::map< Index3, std::vector<BoxKey> >::iterator mi = curdat.fndeidxvec().begin();
                mi != curdat.fndeidxvec().end(); mi++) {
                std::vector<BoxKey>& tmplist = mi->second;
                for (int k = 0; k < tmplist.size(); k++) {
                    BoxKey othkey = tmplist[k];
                    Point3 othctr = BoxCenter(othkey);
                    Point3 tmp = othctr - curctr;
                    tmp /= tmp.l2();
                    Index3 dir = nml2dir(tmp, W);
                    curdat.outdirset().insert(dir);

                    tmp = curctr - othctr;
                    tmp /= tmp.l2();
                    dir = nml2dir(tmp, W);
                    curdat.incdirset().insert(dir);
                }
            }
        }
    }
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::SetupTreeLowFreqLists(BoxKey curkey, BoxDat& curdat) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupTreeLowFreqLists");
#endif
    std::set<BoxKey> Uset, Vset, Wset, Xset;
    // Only deal with boxes strictly in low frequency regime.
    CHECK_TRUE(!IsCellLevelBox(curkey));
    Index3 curpth = curkey.second;
    BoxKey parkey = ParentKey(curkey);
    Index3 parpth = parkey.second;
    //
    int L = pow2(curkey.first);
    int il, iu, jl, ju, kl, ku;
    il = std::max(2 * parpth(0) - 2, 0);  iu = std::min(2 * parpth(0) + 4, L);
    jl = std::max(2 * parpth(1) - 2, 0);  ju = std::min(2 * parpth(1) + 4, L);
    kl = std::max(2 * parpth(2) - 2, 0);  ku = std::min(2 * parpth(2) + 4, L);
    for (int i = il; i < iu; ++i) {
        for (int j = jl; j < ju; ++j) {
            for (int k = kl; k < ku; ++k) {
                Index3 trypth(i, j, k);
                if (trypth(0) == curpth(0)
                    && trypth(1) == curpth(1)
                    && trypth(2) == curpth(2)) {
                    continue;
                }
                BoxKey wntkey(curkey.first, trypth);
                //LEXING: LOOK FOR IT, DO NOT EXIST IF NO CELL BOX COVERING IT
                BoxKey reskey;
                bool found = SetupTreeFind(wntkey, reskey);
                BoxDat& resdat = _boxvec.access(reskey);
                if (!found) {
                    continue;
                }
                bool adj = SetupTreeAdjacent(reskey, curkey);

                if (reskey.first < curkey.first && HasPoints(resdat)) {
                    if (!adj) {
                        Xset.insert(reskey);
                    } else if (IsLeaf(curdat)) {
                        Uset.insert(reskey);
                    }
                }

                if ( reskey.first == curkey.first ) {
                    if (!adj) {
                        Index3 bb = reskey.second - curkey.second;
                        CHECK_TRUE( bb.linfty() <= 3 );
                        if (HasPoints(resdat)) {
                            Vset.insert(reskey);
                        }
                    } else if (IsLeaf(curdat)) {
                        std::queue<BoxKey> rest;
                        rest.push(reskey);
                        while (!rest.empty()) {
                            BoxKey fntkey = rest.front();
                            rest.pop();
                            BoxDat& fntdat = BoxData(fntkey);

                            bool adj = SetupTreeAdjacent(fntkey, curkey);
                            if (!adj && HasPoints(fntdat)) {
                                Wset.insert(fntkey);
                            } 
                            if (adj && IsLeaf(fntdat) && HasPoints(fntdat)) {
                                Uset.insert(fntkey);
                            }
                            if (adj && !IsLeaf(fntdat)) {
                                for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                                    rest.push( ChildKey(fntkey, Index3(CHILD_IND1(ind),
                                                                       CHILD_IND2(ind),
                                                                       CHILD_IND3(ind))) );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (IsLeaf(curdat) && HasPoints(curdat)) {
        Uset.insert(curkey);
    }
    curdat.undeidxvec().insert(curdat.undeidxvec().begin(), Uset.begin(), Uset.end());
    curdat.vndeidxvec().insert(curdat.vndeidxvec().begin(), Vset.begin(), Vset.end());
    curdat.wndeidxvec().insert(curdat.wndeidxvec().begin(), Wset.begin(), Wset.end());
    curdat.xndeidxvec().insert(curdat.xndeidxvec().begin(), Xset.begin(), Xset.end());
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::SetupTreeHighFreqLists(BoxKey curkey, BoxDat& curdat) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupTreeHighFreqLists");
#endif
    Point3 curctr = BoxCenter(curkey);
    double W = BoxWidth(curkey);
    double eps = 1e-12;
    double D = W * W + W; // Far field distance
    double threshold = D - eps;
    if (IsCellLevelBox(curkey)) {
        //LEXING: CHECK THE FOLLOWING
        for (std::map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
             IsCellLevelBox(mi->first); mi++) {
            BoxKey othkey = mi->first;
            BoxDat& othdat = BoxData(othkey);
            if (HasPoints(othdat)) {
                //LEXING: ALWAYS target - source
                Point3 diff = curctr - BoxCenter(othkey);
                if (diff.l2() >= threshold) {
                    Index3 dir = nml2dir(diff / diff.l2(), W);
                    curdat.fndeidxvec()[dir].push_back(othkey);
                } else {
                    curdat.endeidxvec().push_back(othkey);
                }
            }
        }
    } else {
        BoxKey parkey = ParentKey(curkey);
        BoxDat& pardata = BoxData(parkey);
        for (int k = 0; k < pardata.endeidxvec().size(); ++k) {
            BoxKey trykey = pardata.endeidxvec()[k];
            for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                BoxKey othkey = ChildKey(trykey, Index3(CHILD_IND1(ind),
                                                        CHILD_IND2(ind),
                                                        CHILD_IND3(ind)));
                BoxDat& othdat = BoxData(othkey);
                if (HasPoints(othdat)) {
                    //LEXING: ALWAYS target - source
                    Point3 diff = curctr - BoxCenter(othkey);
                    if (diff.l2() >= threshold) {
                        Index3 dir = nml2dir(diff / diff.l2(), W);
                        curdat.fndeidxvec()[dir].push_back(othkey);
                    } else {
                        curdat.endeidxvec().push_back(othkey);
                    }
                }
            }
        }
    }
    return 0;
}

// ----------------------------------------------------------------------
bool Wave3d::SetupTreeFind(BoxKey wntkey, BoxKey& trykey) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupTreeFind");
#endif
    trykey = wntkey;
    while(!IsCellLevelBox(trykey)) {
        std::map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().find(trykey);
        if (mi!=_boxvec.lclmap().end()) {
            return true; //found
        }
        trykey = ParentKey(trykey);
    }
    std::map<BoxKey, BoxDat>::iterator mi = _boxvec.lclmap().find(trykey);
    return (mi != _boxvec.lclmap().end());
}

// ----------------------------------------------------------------------
bool Wave3d::SetupTreeAdjacent(BoxKey meekey, BoxKey youkey) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupTreeAdjacent");
#endif
    int md = std::max(meekey.first,youkey.first);
    Index3 one(1,1,1);
    Index3 meectr(  (2 * meekey.second+one) * pow2(md - meekey.first)  );
    Index3 youctr(  (2 * youkey.second+one) * pow2(md - youkey.first)  );
    int meerad = pow2(md - meekey.first);
    int yourad = pow2(md - youkey.first);
    Index3 dif( ewabs(meectr - youctr) );
    int rad  = meerad + yourad;
    // return true iff at least one edge touch
    return dif[0] <= rad && dif[1] <= rad && dif[2] <= rad && dif.linfty() == rad;
}

//---------------------------------------------------------------------
int Wave3d::setup_Q1(int key, Point3& pos, std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_Q1");
#endif
    int numC = _geomprtn.m();
    Point3 center = ctr();
    Index3 idx;
    for (int d = 0; d < 3; d++) {
        idx(d) = (int) floor(numC * ((pos(d) - center(d) + _K / 2) / _K));
        CHECK_TRUE(idx(d) >= 0 && idx(d) < numC);
    }
    pids.clear();
    pids.push_back( _geomprtn(idx(0),idx(1),idx(2)) ); //JUST ONE PROC
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_Q2(BoxKey boxkey, BoxDat& boxdat, std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_Q2");
#endif
    //for each ent, get all the pids that might need it
    int numC = _geomprtn.m();
    double widC = _K/numC;
    double W = BoxWidth(boxkey);
    if (IsCellLevelBox(boxkey)) {
        //LEXING: CELL LEVEL BOXES ARE NEEDED FOR ALL CPUS
        pids.clear();
        for (int i = 0; i < getMPISize(); i++) {
            pids.push_back(i);
        }
    } else {
        std::set<int> idset;
        Point3 ctr = BoxCenter(boxkey);
        double D = std::max(4 * W * W + 4 * W, 1.0); //LEXING: THIS TAKE CARES THE LOW FREQUENCY PART
        int il = std::max((int)floor((ctr(0)+_K/2-D)/widC),0);
        int iu = std::min((int)ceil( (ctr(0)+_K/2+D)/widC),numC);
        int jl = std::max((int)floor((ctr(1)+_K/2-D)/widC),0);
        int ju = std::min((int)ceil( (ctr(1)+_K/2+D)/widC),numC);
        int kl = std::max((int)floor((ctr(2)+_K/2-D)/widC),0);
        int ku = std::min((int)ceil( (ctr(2)+_K/2+D)/widC),numC);
        //LEXING: IMPROVE THIS
        for (int i = il; i < iu; i++) {
            for (int j = jl; j < ju; j++) {
                for (int k = kl; k < ku; k++) {
                    idset.insert( _geomprtn(i,j,k) );
                }
            }
        }
        pids.clear();
        pids.insert(pids.begin(), idset.begin(), idset.end());
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_Q1_wrapper(int key, Point3& dat, std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_Q1_wrapper");
#endif
    return (Wave3d::_self)->setup_Q1(key, dat, pids);
}

int Wave3d::setup_Q2_wrapper(BoxKey key, BoxDat& dat, std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_Q2_wrapper");
#endif
    return (Wave3d::_self)->setup_Q2(key, dat, pids);
}
