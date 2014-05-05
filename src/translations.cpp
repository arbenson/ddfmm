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

#include <utility>
#include <vector>

int Wave3d::HighFreqM2L(Index3 dir, BoxKey trgkey, DblNumMat& dcp,
			DblNumMat& uep) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2L");
#endif
    Point3 trgctr = BoxCenter(trgkey);
    // get target
    DblNumMat tmpdcp(dcp.m(), dcp.n());
    for (int k = 0; k < tmpdcp.n(); ++k) {
        for (int d = 0; d < 3; ++d) {
            tmpdcp(d, k) = dcp(d, k) + trgctr(d);
        }
    }
    BoxAndDirKey bndkey(trgkey, dir);
    CHECK_TRUE_MSG(_level_prtns.Contains(bndkey, false),
		   "Missing incoming data");
    int mpirank = getMPIRank();
    CHECK_TRUE_MSG(_level_prtns.Owner(bndkey, false) == mpirank,
                   "Updating data that I do not own");
    BoxAndDirDat& bnddat = _level_prtns.Access(bndkey, false);
    CpxNumVec& dcv = bnddat.dirdnchkval();
    for (BoxAndDirKey& key : bnddat.interactionlist()) {
        BoxKey srckey = key._boxkey;
        Point3 srcctr = BoxCenter(srckey);
        // difference vector
        Point3 diff = trgctr - srcctr;
        diff /= diff.l2(); // see wave3d_setup.cpp
	double W = BoxWidth(trgkey);
        CHECK_TRUE( nml2dir(diff, W) == dir );
        // get source
        DblNumMat tmpuep(uep.m(), uep.n());
        for (int k = 0; k < tmpuep.n(); ++k) {
            for (int d = 0; d < 3; ++d) {
                tmpuep(d, k) = uep(d, k) + srcctr(d);
            }
        }
        CHECK_TRUE_MSG(_level_prtns.Contains(key, true), "Missing outgoing data");
        CpxNumVec& ued = _level_prtns.Access(key, true).dirupeqnden();
        CpxNumMat Mts;
        SAFE_FUNC_EVAL( _kernel.kernel(tmpdcp, tmpuep, tmpuep, Mts) );
        // allocate space if necessary
        if (dcv.m() == 0) {
            dcv.resize(Mts.m());
            setvalue(dcv, cpx(0, 0));
        }
        if (ued.m() == 0) {
            ued.resize(Mts.n());
            setvalue(ued, cpx(0, 0));
        }
        SAFE_FUNC_EVAL( zgemv(1.0, Mts, ued, 1.0, dcv) );
    }
    return 0;
}


int Wave3d::HighFreqM2M(BoxAndDirKey& bndkey, NumVec<CpxNumMat>& uc2ue,
                        NumTns<CpxNumMat>& ue2uc) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2M");
#endif
    BoxKey srckey = bndkey._boxkey;
    Index3 dir = bndkey._dir;
    BoxAndDirDat& bnddat = _level_prtns.Access(bndkey, true);
    int mpirank = getMPIRank();
    CHECK_TRUE_MSG(_level_prtns.Owner(bndkey, true) == mpirank,
                   "Updating data that I do not own");
    CpxNumVec& upeqnden = bnddat.dirupeqnden();
    CpxNumVec upchkval(ue2uc(0, 0, 0).m());
    setvalue(upchkval, cpx(0, 0));
    int level = bndkey._boxkey._level;
    CHECK_TRUE_MSG(level <= UnitLevel(), "Incorrect level");

    if (level == UnitLevel()) {
        // The children boxes only have non-directional equivalent densities
        for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
            int a = CHILD_IND1(ind);
            int b = CHILD_IND2(ind);
            int c = CHILD_IND3(ind);
            BoxKey key = ChildKey(srckey, Index3(a, b, c));
            // Do not compute unless we have the child key.
            if (_level_prtns._lf_boxvec.contains(key)) {
                BoxDat& chddat = _level_prtns._lf_boxvec.access(key);
                CHECK_TRUE_MSG(HasPoints(chddat), "No points on child.");
                if (HasPoints(chddat)) {
                    CpxNumVec& chdued = chddat.upeqnden();
                    SAFE_FUNC_EVAL( zgemv(1.0, ue2uc(a, b, c), chdued, 1.0, upchkval) );
                }
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
            BoxAndDirKey bndkey(chdkey, pdir);
            if (_level_prtns.Contains(bndkey, true)) {
                BoxAndDirDat& bnddat = _level_prtns.Access(bndkey, true);
                CpxNumVec& chdued = bnddat.dirupeqnden();
                // TODO(arbenson): fix this
                if (chdued.m() != 0) {
                    SAFE_FUNC_EVAL( zgemv(1.0, ue2uc(a, b, c), chdued, 1.0, upchkval) );
                }
            }
        }
    }

    // Upward check to upward equivalency (uc2ue)
    CpxNumMat& E1 = uc2ue(0);
    CpxNumMat& E2 = uc2ue(1);
    CpxNumMat& E3 = uc2ue(2);
    CpxNumVec tmp0(E3.m());
    CpxNumVec tmp1(E2.m());
    upeqnden.resize(E1.m());
    setvalue(upeqnden, cpx(0, 0));
    SAFE_FUNC_EVAL( zgemv(1.0, E3, upchkval, 0.0, tmp0) );
    SAFE_FUNC_EVAL( zgemv(1.0, E2, tmp0, 0.0, tmp1) );
    SAFE_FUNC_EVAL( zgemv(1.0, E1, tmp1, 0.0, upeqnden) );

    return 0;
}

int Wave3d::HighFreqL2L(Index3 dir, BoxKey trgkey, NumVec<CpxNumMat>& dc2de,
			NumTns<CpxNumMat>& de2dc,
			std::set<BoxAndDirKey>& keys_affected) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqL2L");
#endif
    BoxAndDirKey bndkey(trgkey, dir);
    BoxAndDirDat& bnddat = _level_prtns.Access(bndkey, false);
    CpxNumVec& dnchkval = bnddat.dirdnchkval();
    CpxNumMat& E1 = dc2de(0);
    CpxNumMat& E2 = dc2de(1);
    CpxNumMat& E3 = dc2de(2);
    CpxNumVec tmp0(E3.m());
    CpxNumVec tmp1(E2.m());
    CpxNumVec dneqnden(E1.m());
    
    // TODO(arbenson): remove this check?
    if (dnchkval.m() == 0) {
        return 0;
    }

    CHECK_TRUE_MSG(E3.n() == dnchkval.m(), "E3 mismatch");
    SAFE_FUNC_EVAL( zgemv(1.0, E3, dnchkval, 0.0, tmp0) );
    CHECK_TRUE_MSG(E2.n() == tmp0.m(), "E2 mismatch");
    SAFE_FUNC_EVAL( zgemv(1.0, E2, tmp0, 0.0, tmp1) );
    CHECK_TRUE_MSG(E1.n() == tmp1.m(), "E1 mismatch");
    SAFE_FUNC_EVAL( zgemv(1.0, E1, tmp1, 0.0, dneqnden) );
    dnchkval.resize(0);  // save space

    int level = trgkey._level;
    CHECK_TRUE_MSG(level <= UnitLevel(), "Incorrect level");

    if (level == UnitLevel()) {
        for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
            int a = CHILD_IND1(ind);
            int b = CHILD_IND2(ind);
            int c = CHILD_IND3(ind);
            BoxKey key = ChildKey(trgkey, Index3(a, b, c));
            // If the box was empty, it will not be stored
            if (!_level_prtns._lf_boxvec.contains(key)) {
                continue;
            }
            BoxDat& chddat = _level_prtns._lf_boxvec.access(key);
            CpxNumVec& chddcv = chddat.dnchkval();
            if (chddcv.m() == 0) {
                chddcv.resize(de2dc(a,b,c).m());
                setvalue(chddcv,cpx(0, 0));
            }
            CHECK_TRUE_MSG(de2dc(a, b, c).n() == dneqnden.m(), "Translation mismatch");
            SAFE_FUNC_EVAL( zgemv(1.0, de2dc(a, b, c), dneqnden, 1.0, chddcv) );
        }
    } else {
        Index3 pdir = ParentDir(dir);
        for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
            int a = CHILD_IND1(ind);
            int b = CHILD_IND2(ind);
            int c = CHILD_IND3(ind);             
            BoxKey chdkey = ChildKey(trgkey, Index3(a, b, c));
            BoxAndDirKey bndkey(chdkey, pdir);
            // If we have the data, then we update the directional check values.
            if (_level_prtns.Contains(bndkey, false)) {
                BoxAndDirDat& bnddat = _level_prtns.Access(bndkey, false);
                CpxNumVec& chddcv = bnddat.dirdnchkval();
                if (chddcv.m() == 0) {
                    chddcv.resize(de2dc(a, b, c).m());
                    setvalue(chddcv, cpx(0, 0));
                }
                SAFE_FUNC_EVAL( zgemv(1.0, de2dc(a, b, c), dneqnden, 1.0, chddcv) );
                // We updated the data, so we need to send it back to the children.
                keys_affected.insert(bndkey);
            }
        }
    }
    return 0;
}

int Wave3d::LowFreqM2M(BoxKey& srckey, BoxDat& srcdat, DblNumMat& uep,
                       DblNumMat& ucp, NumVec<CpxNumMat>& uc2ue,
                       NumTns<CpxNumMat>& ue2uc) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqM2M");
#endif
    // TODO(arbenson): What was this being used for before?
    int tdof = 1;

    CHECK_TRUE(HasPoints(srcdat));  // should have points
    Point3 srcctr = BoxCenter(srckey);
    // get array
    CpxNumVec upchkval(tdof * ucp.n());
    setvalue(upchkval, cpx(0, 0));
    CpxNumVec& upeqnden = srcdat.upeqnden();
    // ue2dc
    if (IsLeaf(srcdat)) {
        DblNumMat upchkpos(ucp.m(), ucp.n());
        for (int k = 0; k < ucp.n(); ++k) {
            for (int d = 0; d < dim(); ++d) {
                upchkpos(d, k) = ucp(d, k) + srcctr(d);
            }
        }
        CpxNumMat mat;
        SAFE_FUNC_EVAL( _kernel.kernel(upchkpos, srcdat.extpos(), srcdat.extpos(), mat) );
        SAFE_FUNC_EVAL( zgemv(1.0, mat, srcdat.extden(), 1.0, upchkval) );
    } else {
        for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
            int a = CHILD_IND1(ind);
            int b = CHILD_IND2(ind);
            int c = CHILD_IND3(ind);
            BoxKey key = ChildKey(srckey, Index3(a, b, c));
            if (_level_prtns._lf_boxvec.contains(key)) {
                BoxDat& chddat = _level_prtns._lf_boxvec.access(key);
		CHECK_TRUE(HasPoints(chddat));
                SAFE_FUNC_EVAL( zgemv(1.0, ue2uc(a, b, c), chddat.upeqnden(), 1.0, upchkval) );
            }
        }
    }
    
    // uc2ue
    CpxNumMat& v  = uc2ue(0);
    CpxNumMat& is = uc2ue(1);
    CpxNumMat& up = uc2ue(2);
    CpxNumVec mid(up.m());
    setvalue(mid,cpx(0, 0));
    SAFE_FUNC_EVAL( zgemv(1.0, up, upchkval, 0.0, mid) );
    for (int k = 0; k < mid.m(); ++k) {
        mid(k) = mid(k) * is(k, 0);
    }
    upeqnden.resize(v.m());
    setvalue(upeqnden,cpx(0, 0));
    SAFE_FUNC_EVAL( zgemv(1.0, v, mid, 0.0, upeqnden) );

    return 0;
}

int Wave3d::LowFreqM2L(BoxKey& trgkey, BoxDat& trgdat, DblNumMat& dcp,
                       NumTns<CpxNumTns>& ue2dc, CpxNumVec& dneqnden, DblNumMat& uep,
                       NumVec<CpxNumMat>& dc2de) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqM2L");
#endif
    Point3 trgctr = BoxCenter(trgkey);
    CpxNumVec& dnchkval = trgdat.dnchkval();
    if (dnchkval.m() == 0) {
        dnchkval.resize(dcp.n());
        setvalue(dnchkval,cpx(0, 0));
    }
    if (trgdat.extval().m() == 0) {
        trgdat.extval().resize( trgdat.extpos().n() );
        setvalue(trgdat.extval(), cpx(0, 0));
    }
    DblNumMat dnchkpos(dcp.m(), dcp.n());
    for (int k = 0; k < dcp.n(); ++k) {
        for (int d = 0; d < dim(); ++d) {
            dnchkpos(d, k) = dcp(d, k) + trgctr(d);
        }
    }
    double W = BoxWidth(trgkey);
    // List computations
    SAFE_FUNC_EVAL( UListCompute(trgdat) );
    SAFE_FUNC_EVAL( VListCompute(W, trgdat, trgctr, uep, dcp, dnchkval, ue2dc) );
    SAFE_FUNC_EVAL( WListCompute(W, trgdat, uep) );
    SAFE_FUNC_EVAL( XListCompute(trgdat, dcp, dnchkpos, dnchkval) );
    
    // dnchkval to dneqnden
    CpxNumMat& v  = dc2de(0);
    CpxNumMat& is = dc2de(1);
    CpxNumMat& up = dc2de(2);
    CpxNumVec mid(up.m());
    setvalue(mid,cpx(0, 0));
    SAFE_FUNC_EVAL( zgemv(1.0, up, dnchkval, 0.0, mid) );
    dnchkval.resize(0); // LEXING: SAVE SPACE
    for (int k = 0; k < mid.m(); ++k) {
        mid(k) = mid(k) * is(k, 0);
    }
    dneqnden.resize(v.m());
    setvalue(dneqnden, cpx(0, 0));
    SAFE_FUNC_EVAL( zgemv(1.0, v, mid, 0.0, dneqnden) );
    return 0;
}

int Wave3d::LowFreqL2L(BoxKey& trgkey, BoxDat& trgdat, DblNumMat& dep,
                       NumTns<CpxNumMat>& de2dc, CpxNumVec& dneqnden) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqL2L");
#endif
    Point3 trgctr = BoxCenter(trgkey);

    // Add potentials to children or to exact points
    if (IsLeaf(trgdat)) {
        DblNumMat dneqnpos(dep.m(), dep.n());
        for (int k = 0; k < dep.n(); ++k) {
            for (int d = 0; d < dim(); ++d) {
                dneqnpos(d, k) = dep(d, k) + trgctr(d);
            }
        }
        CpxNumMat mat;
        SAFE_FUNC_EVAL( _kernel.kernel(trgdat.extpos(), dneqnpos, dneqnpos, mat) );
        SAFE_FUNC_EVAL( zgemv(1.0, mat, dneqnden, 1.0, trgdat.extval()) );
    } else {
        // put stuff to children
        for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
            int a = CHILD_IND1(ind);
            int b = CHILD_IND2(ind);
            int c = CHILD_IND3(ind);
            BoxKey key = ChildKey(trgkey, Index3(a, b, c));
            if (!_level_prtns._lf_boxvec.contains(key)) {
                continue;
            }
            BoxDat& chddat = _level_prtns._lf_boxvec.access(key);
            if (chddat.dnchkval().m() == 0) {
                chddat.dnchkval().resize(de2dc(a, b, c).m());
                setvalue(chddat.dnchkval(), cpx(0, 0));
            }
            SAFE_FUNC_EVAL( zgemv(1.0, de2dc(a, b, c), dneqnden, 1.0, chddat.dnchkval()) );
        }
    }

    return 0;
}

int Wave3d::UListCompute(BoxDat& trgdat) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::UListCompute");
#endif
    for (BoxKey& neikey : trgdat.undeidxvec()) {
        BoxDat& neidat = _level_prtns._lf_boxvec.access(neikey);
        CHECK_TRUE(HasPoints(neidat));
        CpxNumMat mat;
        SAFE_FUNC_EVAL( _kernel.kernel(trgdat.extpos(), neidat.extpos(),
				       neidat.extpos(), mat) );
        SAFE_FUNC_EVAL( zgemv(1.0, mat, neidat.extden(), 1.0, trgdat.extval()) );
    }
    return 0;
}

int Wave3d::VListCompute(double W, BoxDat& trgdat, Point3& trgctr, DblNumMat& uep,
                         DblNumMat& dcp, CpxNumVec& dnchkval, NumTns<CpxNumTns>& ue2dc) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::VListCompute");
#endif
    int acc_level = AccLevel();  // log2(1 / \epsilon)
    double step = W / (acc_level - 1);
    setvalue(_valfft, cpx(0, 0));
    for (BoxKey& neikey : trgdat.vndeidxvec()) {
        BoxDat& neidat = _level_prtns._lf_boxvec.access(neikey);
        CHECK_TRUE(HasPoints(neidat));
        Point3 neictr = BoxCenter(neikey);
        Index3 idx;
        for (int d = 0; d < dim(); ++d) {
            idx(d) = int(round( (trgctr[d] - neictr[d]) / W ));
        }
        // Create if it is missing
        if (neidat.fftcnt() == 0) {
            setvalue(_denfft, cpx(0,0));
            CpxNumVec& neiden = neidat.upeqnden();
            for (int k = 0; k < uep.n(); ++k) {
                int a = int( round((uep(0, k) + W / 2) / step) ) + acc_level;
                int b = int( round((uep(1, k) + W / 2) / step) ) + acc_level;
                int c = int( round((uep(2, k) + W / 2) / step) ) + acc_level;
                _denfft(a,b,c) = neiden(k);
            }
            fftw_execute(_fplan);
            neidat.upeqnden_fft() = _denfft;  // COPY to the right place
        }
        CpxNumTns& neidenfft = neidat.upeqnden_fft();
        CpxNumTns& interaction_tensor = ue2dc(idx[0] + 3, idx[1] + 3, idx[2] + 3);
        for (int a = 0; a < 2 * acc_level; ++a) {
            for (int b = 0; b < 2 * acc_level; ++b) {
                for (int c = 0; c < 2 * acc_level; ++c) {
                    _valfft(a, b, c) += (neidenfft(a, b, c) * interaction_tensor(a, b, c));
                }
            }
        }
        // Clean if necessary
        neidat.fftcnt()++;
        if (neidat.fftcnt() == neidat.fftnum()) {
            neidat.upeqnden_fft().resize(0, 0, 0);
            neidat.fftcnt() = 0;
        }
    }
    fftw_execute(_bplan);
    // add back
    double coef = 1.0 / (2 * acc_level * 2 * acc_level * 2 * acc_level);
    for (int k = 0; k < dcp.n(); ++k) {
        int a = int( round((dcp(0, k) + W / 2) / step) ) + acc_level;
        int b = int( round((dcp(1, k) + W / 2) / step) ) + acc_level;
        int c = int( round((dcp(2, k) + W / 2) / step) ) + acc_level;
        dnchkval(k) += (_valfft(a, b, c) * coef);
    }
    return 0;
}

int Wave3d::XListCompute(BoxDat& trgdat, DblNumMat& dcp, DblNumMat& dnchkpos,
                         CpxNumVec& dnchkval) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::XListCompute");
#endif
    for (BoxKey& neikey : trgdat.xndeidxvec()) {
        BoxDat& neidat = _level_prtns._lf_boxvec.access(neikey);
        CHECK_TRUE(HasPoints(neidat));
        Point3 neictr = BoxCenter(neikey);
        if (IsLeaf(trgdat) && trgdat.extpos().n() < dcp.n()) {
            CpxNumMat mat;
            SAFE_FUNC_EVAL( _kernel.kernel(trgdat.extpos(), neidat.extpos(),
					   neidat.extpos(), mat) );
            SAFE_FUNC_EVAL( zgemv(1.0, mat, neidat.extden(), 1.0, trgdat.extval()) );
        } else {
            CpxNumMat mat;
            SAFE_FUNC_EVAL( _kernel.kernel(dnchkpos, neidat.extpos(), neidat.extpos(), mat) );
            SAFE_FUNC_EVAL( zgemv(1.0, mat, neidat.extden(), 1.0, dnchkval) );
        }
    }
    return 0;
}

int Wave3d::WListCompute(double W, BoxDat& trgdat, DblNumMat& uep) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::WListCompute");
#endif
    for (BoxKey& neikey : trgdat.wndeidxvec()) {
        BoxDat& neidat = _level_prtns._lf_boxvec.access(neikey);
        CHECK_TRUE(HasPoints(neidat));
        Point3 neictr = BoxCenter(neikey);
        // upchkpos
        if (IsLeaf(neidat) && neidat.extpos().n() < uep.n()) {
            CpxNumMat mat;
            SAFE_FUNC_EVAL( _kernel.kernel(trgdat.extpos(), neidat.extpos(),
					   neidat.extpos(), mat) );
            SAFE_FUNC_EVAL( zgemv(1.0, mat, neidat.extden(), 1.0, trgdat.extval()) );
        } else {
            double coef = BoxWidth(neikey) / W; // LEXING: SUPER IMPORTANT
            DblNumMat upeqnpos(uep.m(), uep.n()); // local version
            for (int k = 0; k < uep.n(); ++k) {
                for (int d = 0; d < dim(); ++d) {
                    upeqnpos(d, k) = coef * uep(d, k) + neictr(d);
                }
            }
            CpxNumMat mat;
            SAFE_FUNC_EVAL( _kernel.kernel(trgdat.extpos(), upeqnpos, upeqnpos, mat) );
            SAFE_FUNC_EVAL( zgemv(1.0, mat, neidat.upeqnden(), 1.0, trgdat.extval()) );
        }
    }
    return 0;
}
