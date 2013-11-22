/* Distributed Directional Fast Multipole Method
   Copyright (C) 2013 Austin Benson, Lexing Ying, and Jack Poulson

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
#include "mlib3d.hpp"
#include "parallel.hpp"
#include "serialize.hpp"

Point3 ShiftedPoint(int dir_ind, double W) {
#ifndef RELEASE
    CallStackEntry entry("ShiftedPoint");
#endif
    int a = CHILD_IND1(dir_ind);
    int b = CHILD_IND2(dir_ind);
    int c = CHILD_IND3(dir_ind);
    return Point3((a - 0.5) * W / 2, (b - 0.5) * W / 2, (c - 0.5) * W / 2);
}

int Transpose(CpxNumMat& trg, CpxNumMat& src) {
#ifndef RELEASE
    CallStackEntry entry("Transpose");
#endif
    trg.resize(src.n(),src.m());
    for(int i = 0; i < trg.m(); i++) {
        for(int j = 0; j < trg.n(); j++) {
            trg(i, j) = src(j, i);
        }
    }
    return 0;
}

int ApplyShift(DblNumMat& trg, DblNumMat& src, Point3 shift) {
#ifndef RELEASE
    CallStackEntry entry("ApplyShift");
#endif
    for(int k = 0; k < src.n(); k++) {
        for(int d = 0; d < 3; d++) {
            trg(d, k) = src(d, k) + shift(d);
        }
    }
    return 0;
}

int negate(DblNumMat& src) {
#ifndef RELEASE
    CallStackEntry entry("negate");
#endif
    for (int i = 0; i < src.n(); i++) {
        for (int d = 0; d < 3; d++) {
            src(d, i) = -src(d, i);
        }
    }
    return 0;
}

//-----------------------------------
int Mlib3d::setup(std::map<std::string, std::string>& opts)
{
#ifndef RELEASE
    CallStackEntry entry("Mlib3d::setup");
#endif
    //get params
    std::map<std::string, std::string>::iterator mi;
    mi = opts.find("-" + prefix() + "NPQ");
    if(mi!=opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _NPQ;
    }
    mi = opts.find("-" + prefix() + "ldname");
    if(mi!=opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _ldname;
    }
    mi = opts.find("-" + prefix() + "hdname");
    if(mi!=opts.end()) {
        std::istringstream ss(mi->second);
        ss >> _hdname;
    }
  
    //LEXING: read data in a shared way
    std::vector<int> all(1,1);
    std::istringstream liss;
    iC( SharedRead(_ldname, liss) );
    iC( deserialize(_w2ldmap, liss, all) );
    std::istringstream hiss;
    iC( SharedRead(_hdname, hiss) );
    iC( deserialize(_w2hdmap, hiss, all) );
  
    return 0;
}

//-----------------------------------
int Mlib3d::UpwardLowFetch(double W, DblNumMat& uep, DblNumMat& ucp,
                           NumVec<CpxNumMat>& uc2ue, NumTns<CpxNumMat>& ue2uc) {
#ifndef RELEASE
    CallStackEntry entry("Mlib3d::UpwardLowFetch");
#endif
    iA(_w2ldmap.find(W) != _w2ldmap.end());
    LowFreqEntry& le = _w2ldmap[W];
  
    uep = le.uep();
    ucp = le.ucp();
    uc2ue.resize(3);
    for (int i = 0; i < 3; ++i) {
      uc2ue(i) = le.uc2ue()(i);
    }
  
    DblNumMat uepchd = _w2ldmap[W / 2].uep();

    ue2uc.resize(2, 2, 2);
    for (int ind = 0; ind < NUM_CHILDREN; ind++) {
        DblNumMat tmp(uepchd.m(), uepchd.n());
        ApplyShift(tmp, uepchd, ShiftedPoint(ind, W));
        iC( _knl.kernel(ucp, tmp, tmp, ue2uc(CHILD_IND1(ind), CHILD_IND2(ind),
                                             CHILD_IND3(ind))) );
    }
    return 0;
}

//-----------------------------------
int Mlib3d::DownwardLowFetch(double W, DblNumMat& dep, DblNumMat& dcp,
                             NumVec<CpxNumMat>& dc2de, NumTns<CpxNumMat>& de2dc,
                             NumTns<CpxNumTns>& ue2dc, DblNumMat& uep) {
#ifndef RELEASE
    CallStackEntry entry("Mlib3d::DownwardLowFetch");
#endif
    iA(_w2ldmap.find(W) != _w2ldmap.end());
    LowFreqEntry& le = _w2ldmap[W];
  
    dep = le.ucp();
    dcp = le.uep();
    uep = le.uep();
    CpxNumMat tmp0( le.uc2ue()(0) );
    CpxNumMat tmp2( le.uc2ue()(2) );
    dc2de.resize(3);
    Transpose(dc2de(0), tmp2);
    dc2de(1) = le.uc2ue()(1);
    Transpose(dc2de(2), tmp0);
    DblNumMat dcpchd = _w2ldmap[W / 2].uep();
  
    de2dc.resize(2,2,2);
    for (int ind = 0; ind < NUM_CHILDREN; ind++) {
        DblNumMat tmp(dcpchd.m(), dcpchd.n());
        ApplyShift(tmp, dcpchd, ShiftedPoint(ind, W));
        iC( _knl.kernel(tmp, dep, dep, de2dc(CHILD_IND1(ind), CHILD_IND2(ind),
                                             CHILD_IND3(ind))) );
    }
  
    ue2dc.resize(7,7,7);
    for (int a = 0; a < 7; a++) {
        for (int b = 0; b < 7; b++) {
            for (int c = 0; c < 7; c++) {
                if (abs(a - 3) > 1 || abs(b - 3) > 1 || abs(c-3) > 1) {
                    ue2dc(a,b,c) = le.ue2dc()(a,b,c);
                }
            }
        }
    }
    return 0;
}

//-----------------------------------
int Mlib3d::UpwardHighFetch(double W, Index3 dir, DblNumMat& uep, DblNumMat& ucp,
                            NumVec<CpxNumMat>& uc2ue, NumTns<CpxNumMat>& ue2uc) {
#ifndef RELEASE
    CallStackEntry entry("Mlib3d::UpwardHighFetch");
#endif
    iA(_w2hdmap.find(W) != _w2hdmap.end());
    std::map<Index3, HghFreqDirEntry>& curmap = _w2hdmap[W];
  
    Index3 srt, sgn, prm;
    iC( HighFetchIndex3Sort(dir, srt, sgn, prm) );
    iA(curmap.count(srt) != 0);

    HghFreqDirEntry& he = curmap[srt];
    DblNumMat ueptmp( he.uep() );
    DblNumMat ucptmp( he.ucp() );
    iC( HighFetchShuffle(prm, sgn, ueptmp, uep) );
    iC( HighFetchShuffle(prm, sgn, ucptmp, ucp) );
    uc2ue.resize(3);
    for (int i = 0; i < 3; i++) {
        uc2ue(i) = he.uc2ue()(i);
    }
  
    DblNumMat uepchd;
    if (W == 1.0) { //unit box
        iA(_w2ldmap.find(W / 2) != _w2ldmap.end());
        LowFreqEntry& le = _w2ldmap[W / 2];
        uepchd = le.uep();
    } else { //large box
        iA(_w2hdmap.find(W / 2) != _w2hdmap.end());
        std::map<Index3,HghFreqDirEntry>& curmap = _w2hdmap[W / 2];
        
        Index3 pdr = predir(dir);
        Index3 srt, sgn, prm;
        iC( HighFetchIndex3Sort(pdr, srt, sgn, prm) );
        iA(curmap.find(srt) != curmap.end());
        HghFreqDirEntry& he = curmap[srt];
        DblNumMat ueptmp( he.uep() );
        iC( HighFetchShuffle(prm, sgn, ueptmp, uepchd) );
    }
    ue2uc.resize(2,2,2);
    for (int ind = 0; ind < NUM_CHILDREN; ind++) {
        DblNumMat tmp(uepchd.m(), uepchd.n());
        ApplyShift(tmp, uepchd, ShiftedPoint(ind, W));
        iC( _knl.kernel(ucp, tmp, tmp, ue2uc(CHILD_IND1(ind), CHILD_IND2(ind),
                                             CHILD_IND3(ind))) );
    }
  
    return 0;
}

//-----------------------------------
int Mlib3d::DownwardHighFetch(double W, Index3 dir, DblNumMat& dep, DblNumMat& dcp,
                              NumVec<CpxNumMat>& dc2de, NumTns<CpxNumMat>& de2dc,
                              DblNumMat& uep) {
#ifndef RELEASE
    CallStackEntry entry("Mlib3d::DownwardHighFetch");
#endif
    iA(_w2hdmap.find(W) != _w2hdmap.end());
    std::map<Index3,HghFreqDirEntry>& curmap = _w2hdmap[W];
  
    Index3 srt, sgn, prm;
    iC( HighFetchIndex3Sort(dir, srt, sgn, prm) );
    iA(curmap.find(srt) != curmap.end());
    HghFreqDirEntry& he = curmap[srt];
  
    DblNumMat ueptmp( he.uep() );
    DblNumMat ucptmp( he.ucp() );
    iC( HighFetchShuffle(prm, sgn, ucptmp, dep) ); //ucp->dep
    negate(dep);
    iC( HighFetchShuffle(prm, sgn, ueptmp, dcp) ); //uep->dcp
    negate(dcp);
    iC( HighFetchShuffle(prm, sgn, ueptmp, uep) ); //uep->uep
    dc2de.resize(3);
    for (int k = 0; k < 3; k++) {
        CpxNumMat tmp( he.uc2ue()(2 - k) );
        Transpose(dc2de(k), tmp);
    }
    DblNumMat dcpchd;
    if (W == 1.0) { //unit box
        iA(_w2ldmap.find(W / 2) != _w2ldmap.end());
        LowFreqEntry& le = _w2ldmap[W / 2];
        dcpchd = le.uep();
    } else { //large box
        iA(_w2hdmap.find(W / 2) != _w2hdmap.end());
        std::map<Index3,HghFreqDirEntry>& curmap = _w2hdmap[W / 2];
        
        Index3 pdr = predir(dir);
        Index3 srt, sgn, prm;  iC( HighFetchIndex3Sort(pdr, srt, sgn, prm) );
        iA(curmap.find(srt) != curmap.end());
        HghFreqDirEntry& he = curmap[srt];
        DblNumMat ueptmp( he.uep() );
        iC( HighFetchShuffle(prm, sgn, ueptmp, dcpchd) );
        negate(dcpchd);
    }
  
    de2dc.resize(2,2,2);
    for (int ind = 0; ind < NUM_CHILDREN; ind++) {
        DblNumMat tmp(dcpchd.m(), dcpchd.n());
        ApplyShift(tmp, dcpchd, ShiftedPoint(ind, W));
        iC( _knl.kernel(tmp, dep, dep, de2dc(CHILD_IND1(ind), CHILD_IND2(ind),
                                             CHILD_IND3(ind))) );
    }
  
    return 0;
}

//-----------------------------------
Index3 Mlib3d::predir(Index3 dir) {
#ifndef RELEASE
    CallStackEntry entry("Mlib3d::predir");
#endif
    int C = dir.linfty();
    int B = C / 2;
    int midx = -1;
    for (int d = 0; d < 3; d++) {
        if (abs(dir(d)) == C) {
            midx = d;
        }
    }
    //midx gives the direction
    Index3 res;
    for (int d = 0; d < 3; d++) {
        res(d) = (dir(d) + C - 1) / 2;
        res(d) = 2 * (res(d) / 2) + 1 - B;
    }
    res(midx) = dir(midx) / 2;
    return res;
}



//---------------------------------------------------------------------
int Mlib3d::HighFetchShuffle(Index3 prm, Index3 sgn, DblNumMat& tmp,
                             DblNumMat& res) {
#ifndef RELEASE
    CallStackEntry entry("Mlib3d::HighFetchShuffle");
#endif
    res.resize(3, tmp.n());
    for (int k = 0; k < res.n(); k++) {
        for (int d = 0; d < 3; d++) {
            res(prm(d), k) = tmp(d, k);
        }
    }
    for (int k = 0; k < res.n(); k++) {
        for (int d = 0; d < 3; d++) {
            res(d, k) = sgn(d) * res(d, k);
        }
    }
    return 0;
}

//---------------------------------------------------------------------
int Mlib3d::HighFetchIndex3Sort(Index3 val, Index3& srt, Index3& sgn,
                                Index3& prm) {
#ifndef RELEASE
    CallStackEntry entry("Mlib3d::HighFetchIndex3Sort");
#endif
    //make it positive
    for (int d = 0; d < 3; d++) {
        sgn(d) = 1;
        if (val(d) < 0) {
            sgn(d) = -1;
        }
        val(d) = abs(val(d));
    }
    //sort
    int upp = val.linfty() + 1;
    for (int d = 0; d < 3; d++) {
        int minidx = d;
        int minval = val[d];
        for (int e = 0; e < 3; e++) {
            if (val[e] < minval) {
                minval = val[e];
                minidx = e;
            }
        }
        srt[d] = minval;
        prm[d] = minidx;
        val[minidx] = upp;
    }
    return 0;
}


//-------------------
int serialize(const LowFreqEntry& le, std::ostream& os,
              const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    iC( serialize(le._uep, os, mask) );
    iC( serialize(le._ucp, os, mask) );
    iC( serialize(le._uc2ue, os, mask) );
    iC( serialize(le._ue2dc, os, mask) );
    return 0;
}
int deserialize(LowFreqEntry& le, std::istream& is,
                const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    iC( deserialize(le._uep, is, mask) );
    iC( deserialize(le._ucp, is, mask) );
    iC( deserialize(le._uc2ue, is, mask) );
    iC( deserialize(le._ue2dc, is, mask) );
    return 0;
}

//-------------------
int serialize(const HghFreqDirEntry& he, std::ostream& os,
              const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    iC( serialize(he._uep, os, mask) );
    iC( serialize(he._ucp, os, mask) );
    iC( serialize(he._uc2ue, os, mask) );
    return 0;
}

int deserialize(HghFreqDirEntry& he, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    iC( deserialize(he._uep, is, mask) );
    iC( deserialize(he._ucp, is, mask) );
    iC( deserialize(he._uc2ue, is, mask) );
    return 0;
}
