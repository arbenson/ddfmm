#ifndef _MLIB3D_HPP_
#define _MLIB3D_HPP_

#include "comobject.hpp"
#include "vec3t.hpp"
#include "numtns.hpp"
#include "offtns.hpp"
#include "kernel3d.hpp"
#include "vecmatop.hpp"

using std::vector;
using std::pair;
using std::map;

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
  int _NPQ;  //int _P;
  string _ldname;
  string _hdname;
  
  //INTERNAL
  map<double, LowFreqEntry> _w2ldmap;
  map<double, map<Index3, HghFreqDirEntry> > _w2hdmap;
  
public:
  Mlib3d(const string& p);
  ~Mlib3d();
  
  Kernel3d& knl() { return _knl; }
  int& NPQ() { return _NPQ; }
  string& ldname() { return _ldname; }
  string& hdname() { return _hdname; }
  
  map<double, LowFreqEntry>& w2ldmap() { return _w2ldmap; }
  map<double, map<Index3, HghFreqDirEntry> >& w2hdmap() { return _w2hdmap; }
  
  int setup(map<string,string>& opts);
  
  int upward_lowfetch(double W, DblNumMat& uep, DblNumMat& ucp, NumVec<CpxNumMat>& uc2ue,
					  NumTns<CpxNumMat>& ue2uc);
  int dnward_lowfetch(double W, DblNumMat& dep, DblNumMat& dcp, NumVec<CpxNumMat>& dc2de,
					  NumTns<CpxNumMat>& de2dc, NumTns<CpxNumTns>& ue2dc, DblNumMat& uep);
  int upward_hghfetch(double W, Index3 dir, DblNumMat& uep, DblNumMat& ucp, NumVec<CpxNumMat>& uc2ue,
					  NumTns<CpxNumMat>& ue2uc);
  int dnward_hghfetch(double W, Index3 dir, DblNumMat& dep, DblNumMat& dcp, NumVec<CpxNumMat>& dc2de,
					  NumTns<CpxNumMat>& de2dc, DblNumMat& uep);
  
  int hghfetch_shuffle(Index3 prm, Index3 sgn, DblNumMat& tmp, DblNumMat& res);
  int hghfetch_index3sort(Index3 val, Index3& srt, Index3& sgn, Index3& prm);
  //int hghfetch_point3sort(Point3 val, Point3& srt, Index3& sgn, Index3& prm);

  Index3 predir(Index3);
};

//-------------------
int serialize(const LowFreqEntry&, ostream&, const vector<int>&);
int deserialize(LowFreqEntry&, istream&, const vector<int>&);
//-------------------
int serialize(const HghFreqDirEntry&, ostream&, const vector<int>&);
int deserialize(HghFreqDirEntry&, istream&, const vector<int>&);

#endif

