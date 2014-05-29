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
#include "DataCollection.hpp"
#include "wave3d.hpp"

#include <vector>

void Wave3d::ConstructLowFreqMap() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::ConstructLowFreqMap");
#endif
    _ldmap.clear();
    int mpirank = getMPIRank();
    double eps = 1e-12;
    for (auto& kv : _level_prtns._lf_boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        double W = BoxWidth(curkey);
        if (HasPoints(curdat) && _level_prtns.Owner(curkey) == mpirank) {
            CHECK_TRUE(W < 1 - eps);
            _ldmap[W].push_back(curkey);
        }
    }
}

void Wave3d::DeleteEmptyBoxes(std::map<BoxKey, BoxDat>& data) {
    std::vector<BoxKey> to_delete;
    for (auto& kv : data) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        if (!HasPoints(curdat)) {
            to_delete.push_back(curkey);
        }
    }
    for (BoxKey& curkey : to_delete) {
        data.erase(curkey);
    }
}

void Wave3d::Finalize() {
    // Go through low frequency boxes and set things to zero.
    for (auto& kv : _level_prtns._lf_boxvec.lclmap()) {
	BoxDat& dat = kv.second;
	setvalue(dat._extden, cpx(0, 0));
	setvalue(dat._upeqnden, cpx(0, 0));
	setvalue(dat._extval, cpx(0, 0));
	setvalue(dat._dnchkval, cpx(0, 0));
	setvalue(dat._upeqnden_fft, cpx(0, 0));
    }


    // Go through high frequency boxes and set things to zero.
    for (LevelBoxAndDirVec& bndvec : _level_prtns._hf_vecs_out) {
	for (auto& kv : bndvec.lclmap()) {
	    BoxAndDirDat& dat = kv.second;
	    setvalue(dat._dirupeqnden, cpx(0, 0));
	    setvalue(dat._dirdnchkval, cpx(0, 0));
	}
    }

    for (LevelBoxAndDirVec& bndvec : _level_prtns._hf_vecs_inc) {
	for (auto& kv : bndvec.lclmap()) {
	    BoxAndDirDat& dat = kv.second;
	    setvalue(dat._dirupeqnden, cpx(0, 0));
	    setvalue(dat._dirdnchkval, cpx(0, 0));
	}
    }

    for (auto& kv : _level_prtns._unit_vec.lclmap()) {
	BoxAndDirDat& dat = kv.second;
	setvalue(dat._dirupeqnden, cpx(0, 0));
	setvalue(dat._dirdnchkval, cpx(0, 0));
    }
}

int Wave3d::setup(std::map<std::string, std::string>& opts) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup");
#endif
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    int mpirank = getMPIRank();
    // Read optional data
    std::map<std::string, std::string>::iterator mi;
    mi = opts.find("-wave3d_ACCU");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _ACCU;
    }
    mi = opts.find("-wave3d_NPQ");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _NPQ;
    }
    mi = opts.find("-wave3d_K");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _K;
    }
    mi = opts.find("-wave3d_ctr");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        double x, y, z;
        ss >> x >> y >> z;
        _ctr = Point3(x,y,z);
    }
    mi = opts.find("-wave3d_ptsmax");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _ptsmax;
    }
    mi = opts.find("-wave3d_maxlevel");
    if (mi != opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _maxlevel;
    }

    if (mpirank == 0) {
        std::cout << _K <<      " | "
                  << _ACCU <<   " | "
                  << _NPQ <<    " | "
                  << _ctr <<    " | "
                  << _ptsmax << " | "
                  << _maxlevel
                  << std::endl;
    }

    // Create the parvecs
    BoxPrtn bp;
    bp.ownerinfo() = _geomprtn;
    _boxvec.prtn() = bp;
    BoxAndDirPrtn tp;
    tp.ownerinfo() = _geomprtn;
    _bndvec.prtn() = tp;

    // Generate octree
    SAFE_FUNC_EVAL(SetupTree());
    int acc_level = AccLevel();

    // Setup FFT stuff.
    _denfft.resize(2 * acc_level, 2 * acc_level, 2 * acc_level);
    _fplan = fftw_plan_dft_3d(2 * acc_level, 2 * acc_level, 2 * acc_level,
			      (fftw_complex*) (_denfft.data()),
                              (fftw_complex*)(_denfft.data()),
			      FFTW_FORWARD, FFTW_MEASURE);
    CHECK_TRUE(_fplan != NULL);
    setvalue(_denfft, cpx(0, 0));
    _valfft.resize(2 * acc_level, 2 * acc_level, 2 * acc_level);
    _bplan = fftw_plan_dft_3d(2 * acc_level, 2 * acc_level, 2 * acc_level,
                              (fftw_complex*) (_valfft.data()),
                              (fftw_complex*) (_valfft.data()),
			      FFTW_BACKWARD, FFTW_ESTIMATE); 
    CHECK_TRUE(_bplan != NULL);
    setvalue(_valfft, cpx(0, 0));

    double t0 = MPI_Wtime();
    // Delete of empty boxes
    DeleteEmptyBoxes(_boxvec.lclmap());

    // Setup of low and high frequency maps
    _level_prtns.Init(_K);
    GatherLocalKeys();
    PrtnDirections(_level_prtns._hdkeys_out, _level_prtns._hf_vecs_out);
    PrtnDirections(_level_prtns._hdkeys_inc, _level_prtns._hf_vecs_inc);
    PrtnUnitLevel();

    // Gather box data at the unit level for the partitioning of the trees.
    std::vector<int> mask1(BoxDat_Number, 0);
    mask1[BoxDat_tag] = 1;
    mask1[BoxDat_ptidxvec] = 1;
    auto transfer_unit_data = [this] (BoxKey key, BoxDat& dat,
				      std::vector<int>& pids) {
	return TransferUnitLevelData(key, dat, pids);
    };
    SAFE_FUNC_EVAL(_boxvec.getBegin(transfer_unit_data, mask1));
    SAFE_FUNC_EVAL(_boxvec.getEnd(mask1));

    // Now we have the unit level information, so we can setup our part of the tree.
    SetupLowFreqOctree();

    // Remove old boxvec data.
    CleanBoxvec();

    // Delete boxes without points
    DeleteEmptyBoxes(_level_prtns._lf_boxvec.lclmap());
    
    // Gather the interaction lists.
    std::vector<int> mask2(BoxAndDirDat_Number, 0);
    mask2[BoxAndDirDat_interactionlist] = 1;
    auto transfer_data = [this] (BoxAndDirKey key, BoxAndDirDat& dat,
				 std::vector<int>& pids) {
	return TransferBoxAndDirData(key, dat, pids);
    };
    SAFE_FUNC_EVAL(_bndvec.getBegin(transfer_data, mask2));
    SAFE_FUNC_EVAL(_bndvec.getEnd(mask2));
    TransferDataToLevels();
    _bndvec.lclmap().clear();

    // Form data maps needed.
    _level_prtns.FormMaps();
    ConstructLowFreqMap();
    double t1 = MPI_Wtime();
    PrintParData(GatherParData(t0, t1), "Partitioning setup.");
    return 0;
}


int Wave3d::RecursiveBoxInsert(std::queue< std::pair<BoxKey, BoxDat> >& tmpq,
                               bool first_pass) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::RecursiveBoxInsert");
#endif
    double eps = 1e-12;
    while (!tmpq.empty()) {
        std::pair<BoxKey, BoxDat> curent = tmpq.front();
        tmpq.pop();
        BoxKey curkey = curent.first;
        BoxDat& curdat = curent.second;
        //LEXING: VERY IMPORTANT
        if (curdat.ptidxvec().size() > 0) {
            curdat.tag() |=  WAVE3D_PTS;
        }
        // We take an action if we are in the high frequency regime with points OR
        // we are in the low frequency regime with sufficient number of points and
        // not at the max depth.
        // If no action is taken, we are at a leaf node.
        bool action = (curkey._level <= UnitLevel() &&
                       curdat.ptidxvec().size() > 0) ||
                       (static_cast<int>(curdat.ptidxvec().size()) > ptsmax() &&
                        curkey._level < maxlevel() - 1);
        if (action) {
            // Subdivide to get new children
            NumTns<BoxDat> chdboxtns(2, 2, 2);
            Point3 curctr = BoxCenter(curkey);
            for (int tmpidx : curdat.ptidxvec()) {
                Point3 tmp = _positions.access(tmpidx); // get position value
                Index3 idx;
                for (int d = 0; d < 3; ++d) {
                    idx(d) = (tmp(d) >= curctr(d));
                }
                // put points to children
                chdboxtns(idx(0), idx(1), idx(2)).ptidxvec().push_back(tmpidx);
            }
            if ((curkey._level != UnitLevel() && first_pass) || !first_pass) {
                // Put non-empty ones into queue
                for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                    int a = CHILD_IND1(ind);
                    int b = CHILD_IND2(ind);
                    int c = CHILD_IND3(ind);
                    BoxKey key = ChildKey(curkey, Index3(a,b,c));
                    tmpq.push( std::pair<BoxKey, BoxDat>(key, chdboxtns(a, b, c)) );
                }
                // Destory ptidxvec to save memory.  Leave the unit level ones
                // in for later partitioning.
                std::vector<int>().swap(curdat.ptidxvec());
            }
        } else {
            // Copy data into _extpos
            curdat.extpos().resize(3, curdat.ptidxvec().size());
	    if (_kernel.NeedsNormals()) {
		curdat.extnor().resize(3, curdat.ptidxvec().size());
	    }
            for (int g = 0; g < static_cast<int>(curdat.ptidxvec().size()); ++g) {
                int tmpidx = curdat.ptidxvec()[g];
                Point3 tmp = _positions.access(tmpidx);
                for (int d = 0; d < 3; ++d) {
                    curdat.extpos()(d, g) = tmp(d);
                }
		if (_kernel.NeedsNormals()) {
		    // Put in normal vectors.
		    Point3 nor = _normal_vecs.access(tmpidx);
		    for (int d = 0; d < 3; ++d) {
			curdat.extnor()(d, g) = nor(d);
		    }
		}
            }

            //LEXING: VERY IMPORTANT
            curdat.tag() |= WAVE3D_LEAF;
        }
        // Add my self into the tree
        if (first_pass) {
            _boxvec.insert(curkey, curdat);
        } else if (BoxWidth(curkey) < 1 - eps) {
            int mpirank = getMPIRank();
            int owner = _level_prtns.Owner(curkey);
            CHECK_TRUE(mpirank == owner);
            _level_prtns._lf_boxvec.insert(curkey, curdat);
        }
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::SetupTree() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupTree");
#endif
    int mpirank = getMPIRank();
    double K = this->K();

    std::queue< std::pair<BoxKey, BoxDat> > tmpq;
    
    // Get all of the geometry information needed for this processor
    std::vector<int> all(1, 1);
    auto distrib_pts = [this] (int key, Point3& dat, std::vector<int>& pids) {
	return DistribCellPts(key, dat, pids);
    };
    SAFE_FUNC_EVAL(_positions.getBegin(distrib_pts, all));
    SAFE_FUNC_EVAL(_positions.getEnd(all));
    
    int numC = _geomprtn.m();
    int lvlC = CellLevel();
    // Generate cell level boxes, put them into queue
    Point3 bctr = ctr();  // overall center of domain
    NumTns<BoxDat> cellboxtns(numC, numC, numC);
    // Fill boxes with points.
    for (auto& kv : _positions.lclmap()) {
        int key = kv.first;
        Point3 pos = kv.second;
        Index3 idx;
        for (int d = 0; d < 3; ++d) {
             idx(d) = (int) floor(numC * ((pos(d) - bctr(d) + K / 2) / K));
             CHECK_TRUE(idx(d) >= 0 && idx(d) < numC);
        }
        // Put the points in
        cellboxtns(idx(0), idx(1), idx(2)).ptidxvec().push_back( key );
    }
    // Put all boxes owned by this process (whether or not is empty) in a queue.
    // TODO(arbenson): this should be more efficient.
    for (int a = 0; a < numC; ++a) {
        for (int b = 0; b < numC; ++b) {
            for (int c = 0; c < numC; ++c) {
                if (_geomprtn(a, b, c) == mpirank) {
                    BoxKey key(lvlC, Index3(a, b, c));
                    tmpq.push( std::pair<BoxKey, BoxDat>(key, cellboxtns(a, b, c)) );
                }
            }
        }
    }
    cellboxtns.resize(0, 0, 0);

    // Construct the tree.
    RecursiveBoxInsert(tmpq, true);
    SetupHighFreqCallLists();

    // Delete endeidxvec since it was only used to build the interaction lists
    // in the high-frequency regime.
    for (auto& kv : _boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        if (OwnBox(curkey, mpirank)) {
            // Save memory
            std::vector<BoxKey>().swap(curdat.endeidxvec());
        }
    }

    GetHighFreqDirs();
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
    Index3 curpth = curkey._index;
    BoxKey parkey = ParentKey(curkey);
    Index3 parpth = parkey._index;
    //
    int L = pow2(curkey._level);
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
                BoxKey wntkey(curkey._level, trypth);
                // Look for the box.  If it does not exist, no cell box covers it.
                // In this case, we can ignore it.
                BoxKey reskey;
                bool found = SetupTreeFind(wntkey, reskey);
                if (!found) {
                    continue;
                }
                if (!_level_prtns._lf_boxvec.contains(reskey)) {
                    continue;
                }
                BoxDat& resdat = _level_prtns._lf_boxvec.access(reskey);
                bool adj = SetupTreeAdjacent(reskey, curkey);

                if (reskey._level < curkey._level && HasPoints(resdat)) {
                    if (!adj) {
                        Xset.insert(reskey);
                    } else if (IsLeaf(curdat)) {
                        Uset.insert(reskey);
                    }
                }

                if (reskey._level == curkey._level) {
                    if (!adj) {
                        Index3 bb = reskey._index - curkey._index;
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
                            BoxDat& fntdat = _level_prtns._lf_boxvec.access(fntkey);

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
    double D = W * W + W;  // Far field distance
    double threshold = D - eps;
    if (IsCellLevelBox(curkey)) {
        // LEXING: CHECK THE FOLLOWING
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
        for (BoxKey& trykey : pardata.endeidxvec()) {
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
    double eps = 1e-12;
    trykey = wntkey;
    while (BoxWidth(trykey) < 1 - eps) {
        std::map<BoxKey,BoxDat>::iterator mi = _level_prtns._lf_boxvec.lclmap().find(trykey);
        if (mi != _level_prtns._lf_boxvec.lclmap().end()) {
            return true; // found
        }
        trykey = ParentKey(trykey);
    }
    // TODO(arbenson): Fix this.  We need to deal with data in ancestor tree
    // that does not belong to this processor.  For now, it is OK to return always
    // return true, but it is not the most efficient thing to do.
    return true;
#if 0
    std::map<BoxKey, BoxDat>::iterator mi = _level_prtns._lf_boxvec.lclmap().find(trykey);
    return (mi != _level_prtns._lf_boxvec.lclmap().end());
#endif
}

// ----------------------------------------------------------------------
bool Wave3d::SetupTreeAdjacent(BoxKey meekey, BoxKey youkey) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupTreeAdjacent");
#endif
    int md = std::max(meekey._level, youkey._level);
    Index3 one(1, 1, 1);
    Index3 meectr((2 * meekey._index + one) * pow2(md - meekey._level));
    Index3 youctr((2 * youkey._index + one) * pow2(md - youkey._level));
    int meerad = pow2(md - meekey._level);
    int yourad = pow2(md - youkey._level);
    Index3 dif(ewabs(meectr - youctr));
    int rad  = meerad + yourad;
    // return true iff at least one edge touch
    return dif[0] <= rad && dif[1] <= rad && dif[2] <= rad && dif.linfty() == rad;
}

//---------------------------------------------------------------------
int Wave3d::DistribCellPts(int key, Point3& pos, std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::DistribCellPts");
#endif
    int numC = _geomprtn.m();
    Point3 center = ctr();
    Index3 idx;
    for (int d = 0; d < 3; d++) {
        idx(d) = (int) floor(numC * ((pos(d) - center(d) + _K / 2) / _K));
        CHECK_TRUE(idx(d) >= 0 && idx(d) < numC);
    }
    pids.clear();
    // just one processor needs the data
    pids.push_back( _geomprtn(idx(0), idx(1), idx(2)) );
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::DistribUnitPts(int key, Point3& pos, std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::DistribUnitPts");
#endif
    int num_unit_boxes = pow2(UnitLevel());
    Point3 center = ctr();
    Index3 idx;
    for (int d = 0; d < 3; d++) {
        idx(d) = (int) floor(num_unit_boxes * ((pos(d) - center(d) + _K / 2) / _K));
        CHECK_TRUE(idx(d) >= 0 && idx(d) < num_unit_boxes);
    }
    pids.clear();
    int level = UnitLevel();
    Index3 dummy_dir(1, 1, 1);
    BoxAndDirKey bndkey(BoxKey(level, idx), dummy_dir);
    // just one processor needs the data
    pids.push_back( _level_prtns.Owner(bndkey, false) );
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::DistribBoxes(BoxKey boxkey, BoxDat& boxdat, std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::DistribBoxes");
#endif
    int numC = _geomprtn.m();
    double widC = _K / numC;
    double W = BoxWidth(boxkey);
    if (IsCellLevelBox(boxkey)) {
        // Cell level boxes are needed on all processors.
        pids.clear();
        for (int i = 0; i < getMPISize(); ++i) {
            pids.push_back(i);
        }
    } else {
        std::set<int> idset;
        Point3 ctr = BoxCenter(boxkey);
        // This takes care of the low-frequency part.
        double D = std::max(4 * W * W + 4 * W, 1.0);
        int il = std::max((int)floor((ctr(0) + _K / 2 - D) / widC), 0);
        int iu = std::min((int)ceil( (ctr(0) + _K / 2 + D) / widC), numC);
        int jl = std::max((int)floor((ctr(1) + _K / 2 - D) / widC), 0);
        int ju = std::min((int)ceil( (ctr(1) + _K / 2 + D) / widC), numC);
        int kl = std::max((int)floor((ctr(2) + _K / 2 - D) / widC), 0);
        int ku = std::min((int)ceil( (ctr(2) + _K / 2 + D) / widC), numC);
        // TODO(arbenson): improve this
        for (int i = il; i < iu; ++i) {
            for (int j = jl; j < ju; ++j) {
                for (int k = kl; k < ku; ++k) {
                    idset.insert( _geomprtn(i, j, k) );
                }
            }
        }
        pids.clear();
        pids.insert(pids.begin(), idset.begin(), idset.end());
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::DistribLowFreqBoxes(BoxKey boxkey, BoxDat& boxdat, std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::DistribLowFreqBoxes");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    int numC = pow2(UnitLevel());
    double widC = _K / numC;
    double W = BoxWidth(boxkey);
    std::set<int> idset;
    Point3 ctr = BoxCenter(boxkey);
    // This takes care of the low-frequency part.
    double D = std::max(4 * W * W + 4 * W, 1.0);
    int il = std::max((int)floor((ctr(0) + _K / 2 - D) / widC), 0);
    int iu = std::min((int)ceil( (ctr(0) + _K / 2 + D) / widC), numC);
    int jl = std::max((int)floor((ctr(1) + _K / 2 - D) / widC), 0);
    int ju = std::min((int)ceil( (ctr(1) + _K / 2 + D) / widC), numC);
    int kl = std::max((int)floor((ctr(2) + _K / 2 - D) / widC), 0);
    int ku = std::min((int)ceil( (ctr(2) + _K / 2 + D) / widC), numC);
    // TODO(arbenson): improve this
    for (int i = il; i < iu; ++i) {
        for (int j = jl; j < ju; ++j) {
            for (int k = kl; k < ku; ++k) {
                BoxKey key(UnitLevel(), Index3(i, j, k));
                int pid = _level_prtns.Owner(key);
                // The box may not be owned, in which case we just skip it.
                if (pid >= 0 && pid < mpisize) {
                    idset.insert(pid);
                }
            }
        }
    }
    pids.clear();
    pids.insert(pids.begin(), idset.begin(), idset.end());
    return 0;
}

int Wave3d::SetupHighFreqCallLists() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupCallLists");
#endif
    int mpirank = getMPIRank();
    double eps = 1e-12;

    // call get DistribBoxes to get the right boxes
    std::vector<int> mask1(BoxDat_Number, 0);
    mask1[BoxDat_tag] = 1;
    auto distrib_boxes = [this] (BoxKey key, BoxDat& dat,
				 std::vector<int>& pids) {
	return DistribBoxes(key, dat, pids);
    };
    SAFE_FUNC_EVAL(_boxvec.getBegin(distrib_boxes, mask1 ));
    SAFE_FUNC_EVAL(_boxvec.getEnd( mask1 ));
    // Compute lists, low list and high list
    for (auto& kv : _boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        // For all of my boxes with points, setup the call list.
        if (OwnBox(curkey, mpirank) && HasPoints(curdat)) {
            if (BoxWidth(curkey) >= 1 - eps) { // High frequency regime
                SAFE_FUNC_EVAL(SetupTreeHighFreqLists(curkey, curdat));
            }
        }
    }
    return 0;
}

int Wave3d::SetupLowFreqCallLists() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupLowFreqCallLists");
#endif
    double eps = 1e-12;
    int mpirank = getMPIRank();
    std::map<BoxKey, BoxDat>& local_map = _level_prtns._lf_boxvec.lclmap();
    for (auto& kv : local_map) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        if (_level_prtns.Owner(curkey) == mpirank &&
            HasPoints(curdat)) {
            CHECK_TRUE(BoxWidth(curkey) < 1 - eps);
            SAFE_FUNC_EVAL(SetupTreeLowFreqLists(curkey, curdat));
        }
    }
    return 0;
}

int Wave3d::GetExtPos() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::GetExtPos");
#endif
    int mpirank = getMPIRank();
    std::set<BoxKey> reqboxset;
    for (auto& kv : _level_prtns._lf_boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        if (HasPoints(curdat) && _level_prtns.Owner(curkey) == mpirank) {
            reqboxset.insert(curdat.undeidxvec().begin(), curdat.undeidxvec().end());
            reqboxset.insert(curdat.vndeidxvec().begin(), curdat.vndeidxvec().end());
            reqboxset.insert(curdat.wndeidxvec().begin(), curdat.wndeidxvec().end());
            reqboxset.insert(curdat.xndeidxvec().begin(), curdat.xndeidxvec().end());
        }
    }
    std::vector<BoxKey> reqbox;
    reqbox.insert(reqbox.begin(), reqboxset.begin(), reqboxset.end());
    std::vector<int> mask(BoxDat_Number, 0);
    mask[BoxDat_extpos] = 1;
    if (_kernel.NeedsNormals()) {
	mask[BoxDat_extnor] = 1;
    }
    SAFE_FUNC_EVAL( _level_prtns._lf_boxvec.getBegin(reqbox, mask) );
    SAFE_FUNC_EVAL( _level_prtns._lf_boxvec.getEnd(mask) );

    for (auto& kv : _level_prtns._lf_boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        if (HasPoints(curdat) && _level_prtns.Owner(curkey) == mpirank) {
            for (BoxKey& neikey : curdat.vndeidxvec()) {
                BoxDat& neidat = _level_prtns._lf_boxvec.access(neikey);
                neidat.fftnum() ++;
            }
        }
    }
    return 0;
}

int Wave3d::GetHighFreqDirs() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::GetHighFreqDirs");
#endif
    int mpirank = getMPIRank();
    double eps = 1e-12;

    for (auto& kv : _boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
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
            for (auto& kv : curdat.fndeidxvec()) {
                std::vector<BoxKey>& tmplist = kv.second;
                for (int k = 0; k < static_cast<int>(tmplist.size()); ++k) {
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
    return 0;
}


int Wave3d::SetupLowFreqOctree() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::SetupLowFreqOctree");
#endif
    int mpirank = getMPIRank();
    // Put all of the unit level boxes on a queue.
    std::queue< std::pair<BoxKey, BoxDat> > lf_q;
    for (auto& kv : _boxvec.lclmap()) {
        BoxKey key = kv.first;
        BoxDat dat = kv.second;
        int level = key._level;
        Index3 dummy_dir(1, 1, 1);
        BoxAndDirKey bndkey(key, dummy_dir);
        if (level == UnitLevel() &&
            _level_prtns.Owner(bndkey, false) == mpirank) {
            lf_q.push(std::pair<BoxKey, BoxDat>(key, dat));
        }
    }
    
    // Get all of the point information needed.
    std::vector<int> all(1, 1);
    auto distrib_unit_pts = [this] (int key, Point3& dat, std::vector<int>& pids) {
	return DistribUnitPts(key, dat, pids);
    };
    SAFE_FUNC_EVAL( _positions.getBegin(distrib_unit_pts, all) );
    SAFE_FUNC_EVAL( _positions.getEnd(all) );

    // Also get the normal vectors if we need them.
    if (_kernel.NeedsNormals()) {
	std::vector<int> req_nor;
	for (auto& kv : _positions.lclmap()) {
	    req_nor.push_back(kv.first);
	}
	SAFE_FUNC_EVAL( _normal_vecs.getBegin(req_nor, all) );
	SAFE_FUNC_EVAL( _normal_vecs.getEnd(all) );
    }
    
    // Construct the octree in the low-frequency regime.
    RecursiveBoxInsert(lf_q, false);

    std::vector<int> mask(BoxDat_Number, 0);
    mask[BoxDat_tag] = 1;
    auto distrib_low_freq_boxes = [this] (BoxKey key, BoxDat& dat,
					  std::vector<int>& pids) {
	return DistribLowFreqBoxes(key, dat, pids);
    };    
    SAFE_FUNC_EVAL(_level_prtns._lf_boxvec.getBegin(distrib_low_freq_boxes, mask));
    SAFE_FUNC_EVAL(_level_prtns._lf_boxvec.getEnd(mask));

    // Setup the call lists.
    SetupLowFreqCallLists();

    // Get data needed.
    GetExtPos();
    return 0;
}
