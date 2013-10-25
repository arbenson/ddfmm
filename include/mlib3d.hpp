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
#ifndef _MLIB3D_HPP_
#define _MLIB3D_HPP_

#include "comobject.hpp"
#include "vec3t.hpp"
#include "numtns.hpp"
#include "kernel3d.hpp"
#include "vecmatop.hpp"

class LowFreqEntry
{
public:
    DblNumMat _uep;
    DblNumMat _ucp;
    NumVec<CpxNumMat> _uc2ue;
    NumTns<CpxNumTns> _ue2dc;
public:
    LowFreqEntry() {;}
    ~LowFreqEntry() {;}
    DblNumMat& uep() { return _uep; }
    DblNumMat& ucp() { return _ucp; }
    NumVec<CpxNumMat>& uc2ue() { return _uc2ue; }
    NumTns<CpxNumTns>& ue2dc() { return _ue2dc; }
};

class HghFreqDirEntry
{
public:
    DblNumMat _uep;
    DblNumMat _ucp;
    NumVec<CpxNumMat> _uc2ue;
public:
    HghFreqDirEntry() {;}
    ~HghFreqDirEntry() {;}
    DblNumMat& uep() { return _uep; }
    DblNumMat& ucp() { return _ucp; }
    NumVec<CpxNumMat>& uc2ue() { return _uc2ue; }
};

//-----------------------------------
class Mlib3d: public ComObject
{
public:
    //PARAM required
    Kernel3d _knl;
    int _NPQ;
    string _ldname;
    string _hdname;
  
    Mlib3d(const string& p): ComObject(p) {}
    ~Mlib3d() {}
  
    Kernel3d& knl() { return _knl; }
    int& NPQ() { return _NPQ; }
    string& ldname() { return _ldname; }
    string& hdname() { return _hdname; }
  
    map<double, LowFreqEntry>& w2ldmap() { return _w2ldmap; }
    map<double, map<Index3, HghFreqDirEntry> >& w2hdmap() { return _w2hdmap; }
  
    int setup(map<string,string>& opts);
  
    int upwardLowFetch(double W, DblNumMat& uep, DblNumMat& ucp,
                        NumVec<CpxNumMat>& uc2ue, NumTns<CpxNumMat>& ue2uc);
    int downwardLowFetch(double W, DblNumMat& dep, DblNumMat& dcp,
                        NumVec<CpxNumMat>& dc2de, NumTns<CpxNumMat>& de2dc,
                        NumTns<CpxNumTns>& ue2dc, DblNumMat& uep);
    int upwardHighFetch(double W, Index3 dir, DblNumMat& uep, DblNumMat& ucp,
                        NumVec<CpxNumMat>& uc2ue, NumTns<CpxNumMat>& ue2uc);
    int downwardHighFetch(double W, Index3 dir, DblNumMat& dep, DblNumMat& dcp,
                        NumVec<CpxNumMat>& dc2de, NumTns<CpxNumMat>& de2dc,
                        DblNumMat& uep);
  
    int highFetchShuffle(Index3 prm, Index3 sgn, DblNumMat& tmp, DblNumMat& res);
    int highFetchIndex3Sort(Index3 val, Index3& srt, Index3& sgn, Index3& prm);

    Index3 predir(Index3);

private:
    map<double, LowFreqEntry> _w2ldmap;
    map<double, map<Index3, HghFreqDirEntry> > _w2hdmap;
};

//-------------------
int serialize(const LowFreqEntry&, std::ostream&, const std::vector<int>&);
int deserialize(LowFreqEntry&, std::istream&, const std::vector<int>&);
//-------------------
int serialize(const HghFreqDirEntry&, std::ostream&, const std::vector<int>&);
int deserialize(HghFreqDirEntry&, std::istream&, const std::vector<int>&);

#endif
