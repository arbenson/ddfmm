#ifndef _WAVE3D_HPP_
#define _WAVE3D_HPP_

#include "comobject.hpp"
#include "vec3t.hpp"
#include "numtns.hpp"
#include "offtns.hpp"
#include "kernel3d.hpp"
#include "mlib3d.hpp"
#include "parvec.hpp"

using std::vector;
using std::pair;
using std::map;
using std::set;
using std::cerr;
using std::cout;

enum {
  WAVE3D_PTS = 1,
  WAVE3D_TERMINAL = 2,
};

//---------------------------------------------------------------------------

class PtPrtn
{
public:
  vector<int> _ownerinfo;
public:
  PtPrtn() {;}
  ~PtPrtn() {;}
  vector<int>& ownerinfo() { return _ownerinfo; }
  int owner(int key) {
	iA(key<_ownerinfo[_ownerinfo.size()-1]);
	//get the proc which owns the current point
	vector<int>::iterator vi = lower_bound(_ownerinfo.begin(), _ownerinfo.end(), key+1);
	return (vi-_ownerinfo.begin())-1;
  }
};

//---------------------------------------------------------------------------
typedef pair<int,Index3> BoxKey; //level, offset_in_level

class BoxDat
{
public:
  int _tag;
  vector<int> _ptidxvec;
  //
  vector<BoxKey> _undeidxvec;
  vector<BoxKey> _vndeidxvec;
  vector<BoxKey> _wndeidxvec;
  vector<BoxKey> _xndeidxvec;
  vector<BoxKey> _endeidxvec;
  map< Index3, vector<BoxKey> > _fndeidxvec;
  //
  DblNumMat _extpos;
  CpxNumVec _extden;
  CpxNumVec _upeqnden;
  CpxNumVec _extval;
  CpxNumVec _dnchkval;
  //
  CpxNumTns _upeqnden_fft;
  set<Index3> _incdirset;
  set<Index3> _outdirset;
  int _fftnum;
  int _fftcnt;
public:
  BoxDat(): _tag(0), _fftnum(0), _fftcnt(0) {;} //by default, no children
  ~BoxDat() {;}
  //
  int& tag() { return _tag; }
  vector<int>& ptidxvec() { return _ptidxvec; }
  //
  vector<BoxKey>& undeidxvec() { return _undeidxvec; }
  vector<BoxKey>& vndeidxvec() { return _vndeidxvec; }
  vector<BoxKey>& wndeidxvec() { return _wndeidxvec; }
  vector<BoxKey>& xndeidxvec() { return _xndeidxvec; }
  vector<BoxKey>& endeidxvec() { return _endeidxvec; }
  map< Index3, vector<BoxKey> >& fndeidxvec() { return _fndeidxvec; }
  //
  DblNumMat& extpos() { return _extpos; }
  CpxNumVec& extden() { return _extden; }
  CpxNumVec& upeqnden() { return _upeqnden; }
  CpxNumVec& extval() { return _extval; }
  CpxNumVec& dnchkval() { return _dnchkval; }
  //
  CpxNumTns& upeqnden_fft() { return _upeqnden_fft; }
  set<Index3>& incdirset() { return _incdirset; }
  set<Index3>& outdirset() { return _outdirset; }
  int& fftnum() { return _fftnum; }
  int& fftcnt() { return _fftcnt; }
};


#define BoxDat_Number 18
enum {
  BoxDat_tag = 0,
  BoxDat_ptidxvec = 1,
  //
  BoxDat_undeidxvec = 2,
  BoxDat_vndeidxvec = 3,
  BoxDat_wndeidxvec = 4,
  BoxDat_xndeidxvec = 5,
  BoxDat_endeidxvec = 6,
  BoxDat_fndeidxvec = 7,
  //
  BoxDat_extpos = 8,
  BoxDat_extden = 9,
  BoxDat_upeqnden = 10,
  BoxDat_extval = 11,
  BoxDat_dnchkval = 12,
  //
  BoxDat_upeqnden_fft = 13,
  BoxDat_incdirset = 14,
  BoxDat_outdirset = 15,
  BoxDat_fftnum = 16,
  BoxDat_fftcnt = 17,
};

class BoxPrtn{
public:
  IntNumTns _ownerinfo;
public:
  BoxPrtn() {;}
  ~BoxPrtn() {;}
  IntNumTns& ownerinfo() { return _ownerinfo; }
  int owner(BoxKey key) {
    int lvl = key.first;
    Index3 idx = key.second;
    int COEF = pow2(lvl) / _ownerinfo.m();
    idx = idx / COEF;
    return _ownerinfo(idx(0), idx(1), idx(2));
  }
};

//---------------------------------------------------------------------------
typedef pair<BoxKey,Index3> BndKey;

class BndDat{
public:
  CpxNumVec _dirupeqnden;
  CpxNumVec _dirdnchkval;
public:
  BndDat() {;}
  ~BndDat() {;}
  CpxNumVec& dirupeqnden() { return _dirupeqnden; }
  CpxNumVec& dirdnchkval() { return _dirdnchkval; }
};

#define BndDat_Number 2
enum {
  BndDat_dirupeqnden = 0,
  BndDat_dirdnchkval = 1,
};

class BndPrtn{
public:
  IntNumTns _ownerinfo;
public:
  BndPrtn() {;}
  ~BndPrtn() {;}
  IntNumTns& ownerinfo() { return _ownerinfo; }
  int owner(BndKey key) {
    int lvl = key.first.first;
    Index3 idx = key.first.second;
    int COEF = pow2(lvl) / _ownerinfo.m();
    idx = idx / COEF;
    return _ownerinfo(idx(0), idx(1), idx(2));
  }
};

//---------------------------------------------------------------------------
class Wave3d: public ComObject
{
public:
  //-----------------------
  ParVec<int, Point3, PtPrtn>* _posptr;
  Kernel3d _knl;
  int _ACCU;
  int _NPQ;
  Mlib3d* _mlibptr; //read data in parallel and then send to other processors
  IntNumTns _geomprtn;
  //
  double _K;
  Point3 _ctr;
  int _ptsmax;
  int _maxlevel;
  //
  ParVec<BoxKey, BoxDat, BoxPrtn> _boxvec;
  ParVec<BndKey, BndDat, BndPrtn> _bndvec;
  //
  CpxNumTns _denfft, _valfft;
  fftw_plan _fplan, _bplan;
  //
  static Wave3d* _self;
public:
  Wave3d(const string& p);
  ~Wave3d();
  //member access
  ParVec<int, Point3, PtPrtn>*& posptr() { return _posptr; }
  Kernel3d& knl() { return _knl; }
  int& ACCU() { return _ACCU; }
  int& NPQ() { return _NPQ; }
  Mlib3d*& mlibptr() { return _mlibptr; }
  IntNumTns& geomprtn() { return _geomprtn; }
  //
  int P() {	assert(_ACCU>=1 && _ACCU<=3);	if(_ACCU==1) return 4;	else if(_ACCU==2) return 6;	else return 8;  }
  double& K() { return _K; }
  Point3& ctr() { return _ctr; }
  int& ptsmax() { return _ptsmax; }
  int& maxlevel() { return _maxlevel; }
  double width() { return _K; }
  //access information from BoxKey
  Point3 center(BoxKey& curkey) {
    int depth = curkey.first;
    Index3 path = curkey.second;
    int tmp = pow2(depth);
    Point3 t;	for(int d=0; d<3; d++)	  t(d) = _ctr(d) - _K/2 + (path(d)+0.5)/tmp*_K;
    return t;
  }
  double width(BoxKey& curkey) {
    int depth = curkey.first;
    int tmp = pow2(depth);
    return _K/tmp;
  }
  bool iscell(const BoxKey& curkey) {
    return (curkey.first==celllevel());
  }
  BoxKey parkey(BoxKey& curkey) {
    BoxKey tmp;	tmp.first = curkey.first-1;	tmp.second = curkey.second/2;
    return tmp;
  }
  BoxKey chdkey(BoxKey& curkey, Index3 idx) {
    BoxKey tmp;	tmp.first = curkey.first + 1;	tmp.second = 2*curkey.second + idx;
    return tmp;
  }
  BoxDat& boxdata(BoxKey& curkey) {
    return _boxvec.access(curkey);
  }
  bool isterminal(BoxDat& curdat) {
    return (curdat.tag()&WAVE3D_TERMINAL);
  }
  bool ispts(BoxDat& curdat) {
    return (curdat.tag()&WAVE3D_PTS);
  }
  
  //simple functions
  int dim() { return 3; }
  int unitlevel() { return int(round(log(_K)/log(2))); } //the level such that the box has width 1
  int celllevel() { int numC = _geomprtn.m(); return int(round(log(numC)/log(2))); }
  Index3 nml2dir(Point3 nml, double W);
  Index3 predir(Index3 dir);
  vector<Index3> chddir(Index3 dir);
  double dir2width(Index3 dir);
  //main functions
  int setup(map<string,string>& opts);
  int eval( ParVec<int, cpx, PtPrtn>& den, ParVec<int, cpx, PtPrtn>& val);
  int check(ParVec<int, cpx, PtPrtn>& den, ParVec<int, cpx, PtPrtn>& val, IntNumVec& chkkeyvec, double& relerr);
  //aux functions
  int setup_tree();
  static int setup_Q1_wrapper(int key, Point3& dat, vector<int>& pids);
  static int setup_Q2_wrapper(BoxKey key, BoxDat& dat, vector<int>& pids);
  int setup_Q1(int key, Point3& dat, vector<int>& pids);
  int setup_Q2(BoxKey key, BoxDat& dat, vector<int>& pids);
  int setup_tree_callowlist( BoxKey, BoxDat& );
  int setup_tree_calhghlist( BoxKey, BoxDat& );
  bool setup_tree_find(BoxKey wntkey, BoxKey& reskey);
  bool setup_tree_adjacent(BoxKey me, BoxKey yo);
  //
  int eval_upward_low(double W, vector<BoxKey>&, set<BoxKey>& reqboxset); //may be we can make it to be local index
  int eval_dnward_low(double W, vector<BoxKey>&);
  //
  int eval_upward_hgh_recursive(double W, Index3 nowdir, map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap, set<BndKey>& reqbndset);
  int eval_dnward_hgh_recursive(double W, Index3 nowdir, map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap);
  int eval_upward_hgh(double W, Index3 dir, pair< vector<BoxKey>, vector<BoxKey> >& hdvecs, set<BndKey>& reqbndset);
  int eval_dnward_hgh(double W, Index3 dir, pair< vector<BoxKey>, vector<BoxKey> >& hdvecs);
  //
  int mpirank() const { int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); return rank; }
  int mpisize() const { int size; MPI_Comm_size(MPI_COMM_WORLD, &size); return size; }
};

//-------------------
int serialize(const PtPrtn&, ostream&, const vector<int>&);
int deserialize(PtPrtn&, istream&, const vector<int>&);
//-------------------
int serialize(const BoxDat&, ostream&, const vector<int>&);
int deserialize(BoxDat&, istream&, const vector<int>&);
//-------------------
//BoxPrtn, not necessary
//-------------------
int serialize(const BndDat&, ostream&, const vector<int>&);
int deserialize(BndDat&, istream&, const vector<int>&);
//-------------------
//BndPrtn, not necessary


#endif

