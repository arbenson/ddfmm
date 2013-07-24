#ifndef _WAVE3D_HPP_
#define _WAVE3D_HPP_

#include "comobject.hpp"
#include "vec3t.hpp"
#include "numtns.hpp"
#include "kernel3d.hpp"
#include "mlib3d.hpp"
#include "parvec.hpp"

using std::vector;
using std::pair;
using std::map;
using std::set;
using std::cerr;
using std::cout;

#define NUM_DIRS (8)
#define DIR_1(x) ((x & 4) >> 2)
#define DIR_2(x) ((x & 2) >> 1)
#define DIR_3(x) ((x & 1))

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
#ifndef RELEASE
	CallStackEntry entry("PtPrtn::owner");
#endif
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
    // TODO (Austin): Some of these should be private
    int _fftnum;
    int _fftcnt;
    //
    DblNumMat _extpos;  // positions of exact points (leaf level)
    CpxNumVec _extden;  // Exact densities  
    CpxNumVec _upeqnden;  // Upward equivalent density
    CpxNumVec _extval;  // Exact potential value
    CpxNumVec _dnchkval; // Downward check potential

    int _tag;
    vector<int> _ptidxvec;
    //
    vector<BoxKey> _undeidxvec;  // U List
    vector<BoxKey> _vndeidxvec;  // V List
    vector<BoxKey> _wndeidxvec;  // W List
    vector<BoxKey> _xndeidxvec;  // X List
    vector<BoxKey> _endeidxvec;  // Close directions
    map< Index3, vector<BoxKey> > _fndeidxvec; // Far away directions

    // Auxiliarly data structures for FFT
    CpxNumTns _upeqnden_fft;
    set<Index3> _incdirset;
    set<Index3> _outdirset;


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

class BoxPrtn {
public:
    IntNumTns _ownerinfo;
public:
    BoxPrtn() {;}
    ~BoxPrtn() {;}
    IntNumTns& ownerinfo() { return _ownerinfo; }
    int owner(BoxKey key) {
#ifndef RELEASE
	CallStackEntry entry("BoxPrtn::owner");
#endif
        int lvl = key.first;
        Index3 idx = key.second;
        int COEF = pow2(lvl) / _ownerinfo.m();
        idx = idx / COEF;
        return _ownerinfo(idx(0), idx(1), idx(2));
    }
};

//---------------------------------------------------------------------------
typedef pair<BoxKey,Index3> BndKey;

// Boundary data
class BndDat {
public:
    CpxNumVec _dirupeqnden;
    CpxNumVec _dirdnchkval;
public:
    BndDat() {;}
    ~BndDat() {;}
    // Directional upward equivalent density
    CpxNumVec& dirupeqnden() { return _dirupeqnden; }
    // Directional downward check value
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
#ifndef RELEASE
	CallStackEntry entry("BndPrtn::owner");
#endif
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
    Mlib3d* _mlibptr; // Read data in parallel and then send to other processors
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
    double& K() { return _K; }
    Point3& ctr() { return _ctr; }
    int& ptsmax() { return _ptsmax; }
    int& maxlevel() { return _maxlevel; }

    //main functions
    int setup(map<string,string>& opts);

    // Compute the potentials at the target points.
    // den contains the density information
    // TODO (Austin): what is val here?
    int eval( ParVec<int, cpx, PtPrtn>& den, ParVec<int, cpx, PtPrtn>& val);

    // Compute the true solution and store the relative err in relerr.
    int check(ParVec<int, cpx, PtPrtn>& den, ParVec<int, cpx, PtPrtn>& val,
              IntNumVec& chkkeyvec, double& relerr);

private:
    double width() { return _K; }
    //access information from BoxKey
    Point3 center(BoxKey& curkey);
    double width(BoxKey& curkey) { return _K / pow2(curkey.first); }
    bool iscell(const BoxKey& curkey) { return curkey.first == cell_level(); }

    // Return the key of the parent box of the box corresponding to curkey.
    BoxKey parkey(BoxKey& curkey) {
	return BoxKey(curkey.first - 1, curkey.second / 2);
    }

    // Return the key of a child box of the box corresponding to curkey.
    // The index into the 8 children is given by idx
    BoxKey chdkey(BoxKey& curkey, Index3 idx) {
        return BoxKey(curkey.first + 1, 2 * curkey.second + idx);
    }

    BoxDat& boxdata(BoxKey& curkey) { return _boxvec.access(curkey); }

    bool isterminal(BoxDat& curdat) { return curdat.tag() & WAVE3D_TERMINAL; }

    // Returns true iff curdat contains points.
    bool has_pts(BoxDat& curdat) { return curdat.tag() & WAVE3D_PTS; }

    // Determine whether the box corresponding to curkey is owned
    // by this processor.
    bool own_box(BoxKey& curkey, int mpirank) { return _boxvec.prtn().owner(curkey) == mpirank; }

    // Return dimension of this problem.    
    int dim() { return 3; }

    //the level such that the box has width 1
    int unitlevel() { return int(round(log(_K) / log(2))); } 
    int cell_level() { return int(round(log(_geomprtn.m()) / log(2))); }
    Index3 nml2dir(Point3 nml, double W);

    // Compute the parent direction given the child direction.  If the child
    // direction is for a box B, then the parent direction is a direction
    // associated with all children boxes C of B.
    //
    // dir is the child direction
    Index3 predir(Index3 dir);
    vector<Index3> chddir(Index3 dir);
    double dir2width(Index3 dir);

    int setup_tree();
    static int setup_Q1_wrapper(int key, Point3& dat, vector<int>& pids);
    static int setup_Q2_wrapper(BoxKey key, BoxDat& dat, vector<int>& pids);
    int setup_Q1(int key, Point3& dat, vector<int>& pids);
    int setup_Q2(BoxKey key, BoxDat& dat, vector<int>& pids);
    int setup_tree_callowlist( BoxKey, BoxDat& );
    int setup_tree_calhghlist( BoxKey, BoxDat& );
    bool setup_tree_find(BoxKey wntkey, BoxKey& reskey);
    bool setup_tree_adjacent(BoxKey me, BoxKey yo);

    int P();

    pair<double, double> mean_var(time_t t0, time_t t1);


    // Functions for evaluation

    // Travel up the octree and visit the boxes in the low frequency regime.
    // Compute outgoing non-directional equivalent densities using M2M.
    //
    // W is the width of the box
    // srcvec is the list of all boxes owned by this processor of width W
    // reqboxset is filled with all boxes whose information this processor needs
    // for the low-frequency downward pass
    int eval_upward_low(double W, vector<BoxKey>& srcvec, set<BoxKey>& reqboxset);

    // Travel down the octree and visit the boxes in the low frequency regime.
    //
    // W is the width of the box
    int eval_dnward_low(double W, vector<BoxKey>& trgvec);

    //
    int eval_upward_hgh_recursive(double W, Index3 nowdir,
            map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap,
            set<BndKey>& reqbndset);
    int eval_dnward_hgh_recursive(double W, Index3 nowdir,
            map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap);
    
    int eval_upward_hgh(double W, Index3 dir,
            pair< vector<BoxKey>, vector<BoxKey> >& hdvecs, set<BndKey>& reqbndset);
    int eval_dnward_hgh(double W, Index3 dir,
                        pair< vector<BoxKey>, vector<BoxKey> >& hdvecs);
    
    int U_list_compute(BoxDat& trgdat);
    int X_list_compute(BoxDat& trgdat, DblNumMat& dcp, DblNumMat& dnchkpos,
                       CpxNumVec& dnchkval);
    int W_list_compute(BoxDat& trgdat, double W, DblNumMat& uep);
    int V_list_compute(BoxDat& trgdat, double W, int _P, Point3& trgctr,
                       DblNumMat& uep, DblNumMat& dcp, CpxNumVec& dnchkval,
                       NumTns<CpxNumTns>& ue2dc);

    int get_reqs(Index3 dir, pair< vector<BoxKey>, vector<BoxKey> >& hdvecs,
		 set<BndKey>& reqbndset);
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
