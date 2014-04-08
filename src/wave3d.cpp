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

//---------------------------------------------------------------------------
Wave3d* Wave3d::_self = NULL;

//-----------------------------------
Wave3d::Wave3d(const std::string& p): ComObject(p), _posptr(NULL), _mlibptr(NULL),
                                      _fplan(NULL), _bplan(NULL), _ACCU(1), _NPQ(4),
			              _K(64), _ctr(Point3(0, 0, 0)), _ptsmax(100) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::Wave3d");
#endif
    _self = this;
}
//-----------------------------------
Wave3d::~Wave3d() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::~Wave3d");
#endif
    _self = this;
    if ( _fplan != NULL) { fftw_destroy_plan(_fplan); }
    if ( _bplan != NULL) { fftw_destroy_plan(_bplan); }
}

//-----------------------------------
Index3 Wave3d::nml2dir(Point3 n, double W) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::nml2dir");
#endif
    int C = _NPQ * int(round(W));
    int B = C / 2;
    int midx = 0;
    double mval = abs(n(0));
    for (int d = 1; d < 3; ++d) {
	if (mval < abs(n(d))) {
            midx = d;
            mval = abs(n(d));
	}
    }
    //midx gives the direction (can be + or -)
    Point3 val = n / mval;
    for (int d = 0; d < 3; ++d) {
        val(d) = atan(val(d));
    }
    double Ang = (M_PI/2)/C;
    Index3 res;
    for (int d = 0; d < 3; ++d) {
        int tmp = int(floor(val(d) / Ang));
        tmp = std::min(std::max(tmp, -B), B - 1);
        res(d) = 2 * tmp + 1;
    }
    res(midx) = C * int(round(val(midx))); //val(midx)==1 or -1
    return res;
}

//-----------------------------------
Index3 Wave3d::ParentDir(Index3 dir) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::ParentDir");
#endif
    int C = dir.linfty();
    int B = C / 2;
    int midx = -1;
    for (int d = 0; d < 3; ++d) {
      if (abs(dir(d)) == C) {
	midx = d;
      }
    }
    assert(midx!=-1);
    //midx gives the direction
    Index3 res;
    for (int d = 0; d < 3; ++d) {
        res(d) = (dir(d) + C - 1) / 2;
        res(d) = 2 * (res(d) / 2) + 1 - B;
    }
    res(midx) = dir(midx) / 2;
    return res;
}

//-----------------------------------
std::vector<Index3> Wave3d::ChildDir(Index3 dir) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::ChildDir");
#endif
    int C = dir.linfty();
    std::vector<int> oidx;
    for (int d = 0; d < 3; ++d) {
	if(abs(dir(d)) != C) {
	    oidx.push_back(d);
	}
    }
    std::vector<Index3> res;
    for (int a = 0; a < 2; ++a) {
	for (int b = 0; b < 2; ++b) {
	    Index3 tmp = 2 * dir;
	    tmp(oidx[0]) += 2 * a - 1;
	    tmp(oidx[1]) += 2 * b - 1;
	    res.push_back(tmp);
	}
    }
    return res;
}

//-----------------------------------
double Wave3d::Dir2Width(Index3 dir) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::dir2width");
#endif
    int C = dir.linfty();
    return double(C / _NPQ);
}

//-----------------------------------------------------------
Point3 Wave3d::BoxCenter(BoxKey& curkey) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::center");
#endif
    Index3 path = curkey.second;
    int tmp = pow2(curkey.first);
    Point3 t;
    for (int d = 0; d < 3; ++d) {
	t(d) = _ctr(d) - _K / 2 + (path(d) + 0.5) / tmp * _K;
    }
    return t;
}

//-----------------------------------------------------------
int Wave3d::P() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::P");
#endif
    CHECK_TRUE(_ACCU >= 1 && _ACCU <= 3);
    switch (_ACCU) {
    case 1:
	return 4;
    case 2:
	return 6;
    default:
	return 8;
    }
}

//--------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------
int serialize(const PtPrtn& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    serialize(val._ownerinfo, os, mask);
    return 0;
}

//-----------------------------------------------------------
int deserialize(PtPrtn& val, std::istream& is, const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    deserialize(val._ownerinfo, is, mask);
    return 0;
}

//-----------------------------------------------------------
int serialize(const BoxDat& val, std::ostream& os, const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int i = 0;
  
    if (mask[i] == 1) serialize(val._tag, os, mask);  i++;
    if (mask[i] == 1) serialize(val._ptidxvec, os, mask);  i++;
  
    if (mask[i] == 1) serialize(val._undeidxvec, os, mask);  i++;
    if (mask[i] == 1) serialize(val._vndeidxvec, os, mask);  i++;
    if (mask[i] == 1) serialize(val._wndeidxvec, os, mask);  i++;
    if (mask[i] == 1) serialize(val._xndeidxvec, os, mask);  i++;
    if (mask[i] == 1) serialize(val._endeidxvec, os, mask);  i++;
    if (mask[i] == 1) serialize(val._fndeidxvec, os, mask);  i++;
  
    if (mask[i] == 1) serialize(val._extpos, os, mask);  i++;
    if (mask[i] == 1) serialize(val._extden, os, mask);  i++;
    if (mask[i] == 1) serialize(val._upeqnden, os, mask);  i++;
    if (mask[i] == 1) serialize(val._extval, os, mask);  i++;
    if (mask[i] == 1) serialize(val._dnchkval, os, mask);  i++;
  
    if (mask[i] == 1) serialize(val._upeqnden_fft, os, mask);  i++;
    if (mask[i] == 1) serialize(val._incdirset, os, mask);  i++;
    if (mask[i] == 1) serialize(val._outdirset, os, mask);  i++;
    if (mask[i] == 1) serialize(val._fftnum, os, mask);  i++;
    if (mask[i] == 1) serialize(val._fftcnt, os, mask);  i++;
  
    CHECK_TRUE(i == BoxDat_Number);
  
    return 0;
}

//-----------------------------------------------------------
int deserialize(BoxDat& val, std::istream& is, const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int i = 0;
  
    if (mask[i] == 1) deserialize(val._tag, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._ptidxvec, is, mask);  i++;

    if (mask[i] == 1) deserialize(val._undeidxvec, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._vndeidxvec, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._wndeidxvec, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._xndeidxvec, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._endeidxvec, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._fndeidxvec, is, mask);  i++;

    if (mask[i] == 1) deserialize(val._extpos, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._extden, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._upeqnden, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._extval, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._dnchkval, is, mask);  i++;
  
    if (mask[i] == 1) deserialize(val._upeqnden_fft, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._incdirset, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._outdirset, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._fftnum, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._fftcnt, is, mask);  i++;
  
    CHECK_TRUE(i == BoxDat_Number);
  
    return 0;
}

//-----------------------------------------------------------
int serialize(const BoxAndDirDat& val, std::ostream& os,
              const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int i = 0;
    if (mask[i] == 1) serialize(val._dirupeqnden, os, mask);  i++;
    if (mask[i] == 1) serialize(val._dirdnchkval, os, mask);  i++;
    CHECK_TRUE(i == BoxAndDirDat_Number);
    return 0;
}
//-----------------------------------------------------------
int deserialize(BoxAndDirDat& val, std::istream& is,
                const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int i = 0;
    if (mask[i] == 1) deserialize(val._dirupeqnden, is, mask);  i++;
    if (mask[i] == 1) deserialize(val._dirdnchkval, is, mask);  i++;
    CHECK_TRUE(i == BoxAndDirDat_Number);
    return 0;
}

//-----------------------------------------------------------
int serialize(const BoxAndDirKey& key, std::ostream& os, const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int i = 0;
    serialize(key._boxkey, os, mask); i++;
    serialize(key._dir, os, mask); i++;
    CHECK_TRUE(i == BoxAndDirKey_Number);
    return 0;
}
//-----------------------------------------------------------
int deserialize(BoxAndDirKey& key, std::istream& is,
                const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int i = 0;
    deserialize(key._boxkey, is, mask);  i++;
    deserialize(key._dir, is, mask);  i++;
    CHECK_TRUE(i == BoxAndDirKey_Number);
    return 0;
}
