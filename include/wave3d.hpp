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
#ifndef _WAVE3D_HPP_
#define _WAVE3D_HPP_

#include "comobject.hpp"
#include "kernel3d.hpp"
#include "mlib3d.hpp"
#include "numtns.hpp"
#include "parvec.hpp"
#include "vec3t.hpp"

#include <algorithm>
#include <vector>

#define NUM_CHILDREN (8)
#define CHILD_IND1(x) ((x & 4) >> 2)
#define CHILD_IND2(x) ((x & 2) >> 1)
#define CHILD_IND3(x) ((x & 1))

enum {
    WAVE3D_PTS = 1,
    WAVE3D_LEAF = 2,
};

//---------------------------------------------------------------------------
class PtPrtn {
public:
    std::vector<int> _ownerinfo;
public:
    PtPrtn() {;}
    ~PtPrtn() {;}
    std::vector<int>& ownerinfo() { return _ownerinfo; }
    int owner(int key) {
#ifndef RELEASE
        CallStackEntry entry("PtPrtn::owner");
#endif
        CHECK_TRUE(key < _ownerinfo[_ownerinfo.size() - 1]);
        // Get the process which owns the current point
        std::vector<int>::iterator vi = lower_bound(_ownerinfo.begin(),
                                                    _ownerinfo.end(), key + 1);
        return (vi - _ownerinfo.begin()) - 1;
    }
};

// Returns true if and only if lhs > rhs in the z-order curve.
bool MortonOrderGreater(Index3 lhs, Index3 rhs);

//---------------------------------------------------------------------------
class BoxKey {
public:
    BoxKey() {;}
    BoxKey(int level, Index3 index) : _level(level), _index(index) {;}
    BoxKey(const BoxKey& key) : _level(key._level), _index(key._index) {;}
    ~BoxKey() {;}

    int _level;
    Index3 _index;

    inline bool operator==(const BoxKey& rhs) const {
        return _level == rhs._level && _index == rhs._index;
    }

    inline bool operator!=(const BoxKey& rhs) const {
        return !operator==(rhs);
    }

    inline bool operator>(const BoxKey& rhs) const {
        return _level > rhs._level ||
              (_level == rhs._level && MortonOrderGreater(_index, rhs._index));
    }

    inline bool operator<(const BoxKey& rhs) const {
        return _level < rhs._level ||
               (_level == rhs._level && MortonOrderGreater(rhs._index, _index));
    }

    inline bool operator>=(const BoxKey& rhs) const {
        return operator>(rhs) || operator==(rhs);
    }

    inline bool operator<=(const BoxKey& rhs) const {
        return operator<(rhs) || operator==(rhs);
    }
};

#define BoxKey_Number 2
enum {
    BoxKey_level = 0,
    BoxKey_index = 1,
};


class BoxDat {
public:
    BoxDat(): _fftnum(0), _fftcnt(0), _tag(0) {;}  // by default, no children
    ~BoxDat() {;}

    // TODO (Austin): Some of these should be private
    int _fftnum;
    int _fftcnt;

    DblNumMat _extpos;   // positions of exact points (leaf level)
    CpxNumVec _extden;   // Exact densities  
    CpxNumVec _upeqnden; // Upward equivalent density
    CpxNumVec _extval;   // Exact potential value
    CpxNumVec _dnchkval; // Downward check potential

    int _tag;
    std::vector<int> _ptidxvec;

    std::vector<BoxKey> _undeidxvec;  // U List
    std::vector<BoxKey> _vndeidxvec;  // V List
    std::vector<BoxKey> _wndeidxvec;  // W List
    std::vector<BoxKey> _xndeidxvec;  // X List
    std::vector<BoxKey> _endeidxvec;  // Boxes in near field
    // _fndeidxvec maps a direction to a vector of boxes that are in the
    // interaction list of this box in that direction
    std::map< Index3, std::vector<BoxKey> > _fndeidxvec;

    // Auxiliarly data structures for FFT
    CpxNumTns _upeqnden_fft;
    std::set<Index3> _incdirset;
    std::set<Index3> _outdirset;

    int& tag() { return _tag; }
    std::vector<int>& ptidxvec() { return _ptidxvec; }

    std::vector<BoxKey>& undeidxvec() { return _undeidxvec; }
    std::vector<BoxKey>& vndeidxvec() { return _vndeidxvec; }
    std::vector<BoxKey>& wndeidxvec() { return _wndeidxvec; }
    std::vector<BoxKey>& xndeidxvec() { return _xndeidxvec; }
    std::vector<BoxKey>& endeidxvec() { return _endeidxvec; }
    std::map< Index3, std::vector<BoxKey> >& fndeidxvec() { return _fndeidxvec; }

    DblNumMat& extpos() { return _extpos; }
    CpxNumVec& extden() { return _extden; }
    CpxNumVec& upeqnden() { return _upeqnden; }
    CpxNumVec& extval() { return _extval; }
    CpxNumVec& dnchkval() { return _dnchkval; }

    CpxNumTns& upeqnden_fft() { return _upeqnden_fft; }
    std::set<Index3>& incdirset() { return _incdirset; }
    std::set<Index3>& outdirset() { return _outdirset; }
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
        int lvl = key._level;
        Index3 idx = key._index;
        int COEF = pow2(lvl) / _ownerinfo.m();
        idx = idx / COEF;
        return _ownerinfo(idx(0), idx(1), idx(2));
    }
};

class BoxAndDirKey {
public:
    BoxAndDirKey(BoxKey boxkey, Index3 dir) : _boxkey(boxkey), _dir(dir) {}
    BoxAndDirKey() {}
    ~BoxAndDirKey() {}

    BoxKey _boxkey;
    Index3 _dir;

    inline bool operator==(const BoxAndDirKey& rhs) const {
        return _boxkey == rhs._boxkey && _dir == rhs._dir;
    }

    inline bool operator!=(const BoxAndDirKey& rhs) const {
        return !operator==(rhs);
    }

    inline bool operator>(const BoxAndDirKey& rhs) const {
        return _boxkey > rhs._boxkey ||
               (_boxkey == rhs._boxkey && MortonOrderGreater(_dir, rhs._dir));
    }

    inline bool operator<(const BoxAndDirKey& rhs) const {
        return _boxkey < rhs._boxkey ||
               (_boxkey == rhs._boxkey && MortonOrderGreater(rhs._dir, _dir));
    }

    inline bool operator>=(const BoxAndDirKey& rhs) const {
        return operator>(rhs) || operator==(rhs);
    }

    inline bool operator<=(const BoxAndDirKey& rhs) const {
        return operator<(rhs) || operator==(rhs);
    }
};

#define BoxAndDirKey_Number 2
enum {
    BoxAndDirKey_boxkey = 0,
    BoxAndDirKey_dir = 1,
};


//---------------------------------------------------------------------------
class BoxAndDirLevelPrtn {
public:
    BoxAndDirLevelPrtn() {}
    ~BoxAndDirLevelPrtn() {}
    
    std::vector<BoxAndDirKey> partition_;
    std::vector<BoxAndDirKey> end_partition_;  // for debugging

    // Return process that owns the key.
    int owner(BoxAndDirKey& key);
};

class UnitLevelBoxPrtn {
public:
    UnitLevelBoxPrtn() {}
    ~UnitLevelBoxPrtn() {}
    std::vector<BoxKey> partition_;
    std::vector<BoxKey> end_partition_;  // for debugging

    // Return process that owns the key.
    int owner(BoxAndDirKey& key);
};

class LowFreqBoxPrtn {
public:
    LowFreqBoxPrtn() {;}
    ~LowFreqBoxPrtn() {;}
    std::vector<BoxKey> partition_;
    std::vector<BoxKey> end_partition_;  // for debugging
    int unit_level_;

    // Return process that owns the key.
    int owner(BoxKey& key);
};


//---------------------------------------------------------------------------
// Boundary data
class BoxAndDirDat {
public:
    CpxNumVec _dirupeqnden;
    CpxNumVec _dirdnchkval;
    std::vector<BoxAndDirKey> _interactionlist;

    BoxAndDirDat() {;}
    ~BoxAndDirDat() {;}
    // Directional upward equivalent density
    CpxNumVec& dirupeqnden() { return _dirupeqnden; }
    // Directional downward check value
    CpxNumVec& dirdnchkval() { return _dirdnchkval; }
    // High-frequency interaction list
    std::vector<BoxAndDirKey>& interactionlist() { return _interactionlist; }
};

#define BoxAndDirDat_Number 3
enum {
    BoxAndDirDat_dirupeqnden = 0,
    BoxAndDirDat_dirdnchkval = 1,
    BoxAndDirDat_interactionlist = 2,
};

class BoxAndDirPrtn {
public:
    IntNumTns _ownerinfo;
public:
    BoxAndDirPrtn() {;}
    ~BoxAndDirPrtn() {;}
    IntNumTns& ownerinfo() { return _ownerinfo; }
    int owner(BoxAndDirKey key) {
#ifndef RELEASE
        CallStackEntry entry("BoxAndDirPrtn::owner");
#endif
        int lvl = key._boxkey._level;
        Index3 idx = key._boxkey._index;
        int COEF = pow2(lvl) / _ownerinfo.m();
        idx = idx / COEF;
        return _ownerinfo(idx(0), idx(1), idx(2));
    }
};

typedef ParVec<BoxAndDirKey, BoxAndDirDat, BoxAndDirLevelPrtn> LevelBoxAndDirVec;
typedef std::map< double, std::vector<BoxKey> > ldmap_t;
typedef std::vector< std::vector<BoxAndDirKey> > level_hdkeys_t;
typedef std::vector< std::map<Index3, std::vector<BoxKey> > > level_hdkeys_map_t;

//---------------------------------------------------------------------------
// Handler of all partitions
class LevelPartitions {
public:
  LevelPartitions() {;}
  ~LevelPartitions() {;}

  // out specifies outgoing / incoming
  BoxAndDirDat& Access(BoxAndDirKey key, bool out);
  bool Contains(BoxAndDirKey key, bool out);
  int Owner(BoxAndDirKey key, bool out);
  int Owner(BoxKey key);
  void Init(int K);
  void FormMaps();
  
  std::vector<LevelBoxAndDirVec> _hf_vecs_out;  // outgoing partition for M2M
  std::vector<LevelBoxAndDirVec> _hf_vecs_inc;  // outgoing partition for M2L + L2L

  ParVec<BoxAndDirKey, BoxAndDirDat, UnitLevelBoxPrtn> _unit_vec;
  ParVec<BoxKey, BoxDat, LowFreqBoxPrtn> _lf_boxvec;  // boxes in low frequency regime

  level_hdkeys_t _hdkeys_out;  // which keys I am responsible for (outgoing)
  level_hdkeys_t _hdkeys_inc;  // which keys I am responsible for (incoming)

  level_hdkeys_map_t _level_hdmap_out;
  level_hdkeys_map_t _level_hdmap_inc;

private:
  int unit_level_;
};


//---------------------------------------------------------------------------
class Wave3d: public ComObject {
public:
    ParVec<int, Point3, PtPrtn>* _posptr;
    Kernel3d _kernel;
    int _ACCU;
    int _NPQ;
    Mlib3d* _mlibptr;  // Read data in parallel and then send to other processors
    IntNumTns _geomprtn;

    double _K;
    Point3 _ctr;
    int _ptsmax;
    int _maxlevel;

    ParVec<BoxKey, BoxDat, BoxPrtn> _boxvec;
    ParVec<BoxAndDirKey, BoxAndDirDat, BoxAndDirPrtn> _bndvec;

    CpxNumTns _denfft, _valfft;
    fftw_plan _fplan, _bplan;

    static Wave3d* _self;

    Wave3d(const std::string& p);
    ~Wave3d();
    ParVec<int, Point3, PtPrtn>*& posptr() { return _posptr; }
    Kernel3d& kernel() { return _kernel; }
    int& ACCU() { return _ACCU; }
    int& NPQ() { return _NPQ; }
    Mlib3d*& mlibptr() { return _mlibptr; }
    IntNumTns& geomprtn() { return _geomprtn; }
    double& K() { return _K; }
    Point3& ctr() { return _ctr; }
    int& ptsmax() { return _ptsmax; }
    int& maxlevel() { return _maxlevel; }

    int setup(std::map<std::string, std::string>& opts);

    bool HasPoints(BoxDat& dat) const { return dat.tag() & WAVE3D_PTS; }
    bool IsLeaf(BoxDat& dat) const { return dat.tag() & WAVE3D_LEAF; }

    // Compute the potentials at the target points.
    int eval( ParVec<int, cpx, PtPrtn>& den, ParVec<int, cpx, PtPrtn>& val);

    // Compute the true solution at a few points and return the relative error.
    double check(ParVec<int, cpx, PtPrtn>& den, ParVec<int, cpx, PtPrtn>& val,
                 IntNumVec& chkkeyvec);

private:
    LevelPartitions _level_prtns;
    int _starting_level;

    double width() { return _K; }
    //access information from BoxKey
    Point3 BoxCenter(BoxKey& curkey);
    double BoxWidth(BoxKey& curkey) { return _K / pow2(curkey._level); }
    bool IsCellLevelBox(const BoxKey& curkey) { return curkey._level == CellLevel(); }

    // Return the key of the parent box of the box corresponding to curkey.
    BoxKey ParentKey(BoxKey& curkey) {
        return BoxKey(curkey._level - 1, curkey._index / 2);
    }

    // Return the key of a child box of the box corresponding to curkey.
    // The index into the 8 children is given by idx
    BoxKey ChildKey(BoxKey& curkey, Index3 idx) {
        return BoxKey(curkey._level + 1, 2 * curkey._index + idx);
    }

    BoxDat& BoxData(BoxKey& curkey) { return _boxvec.access(curkey); }

    // Determine whether the box corresponding to curkey is owned
    // by this processor.
    bool OwnBox(BoxKey& curkey, int mpirank) {
        return _boxvec.prtn().owner(curkey) == mpirank;
    }

    // Return dimension of this problem.    
    int dim() { return 3; }

    // The level such that the box has width 1
    int UnitLevel() { return static_cast<int>(round(log(_K) / log(2))); } 

    // The level where the geometry is partitioned.
    int CellLevel() { return static_cast<int>(round(log(_geomprtn.m()) / log(2))); }

    Index3 nml2dir(Point3 nml, double W);

    // Compute the parent direction given the child direction.  If the child
    // direction is for a box B, then the parent direction is a direction
    // associated with all children boxes C of B.
    //
    // dir is the child direction
    Index3 ParentDir(Index3 dir);
    std::vector<Index3> ChildDir(Index3 dir);
    double Dir2Width(Index3 dir);

    // Stuff for setup
    int SetupTree();
    static int DistribCellPts_wrapper(int key, Point3& dat, std::vector<int>& pids);
    static int DistribBoxes_wrapper(BoxKey key, BoxDat& dat, std::vector<int>& pids);
    static int DistribLowFreqBoxes_wrapper(BoxKey key, BoxDat& dat, std::vector<int>& pids);
    static int DistribUnitPts_wrapper(int key, Point3& dat, std::vector<int>& pids);
    int DistribCellPts(int key, Point3& dat, std::vector<int>& pids);
    int DistribBoxes(BoxKey key, BoxDat& dat, std::vector<int>& pids);
    int DistribLowFreqBoxes(BoxKey key, BoxDat& dat, std::vector<int>& pids);
    int DistribUnitPts(int key, Point3& dat, std::vector<int>& pids);
    int SetupTreeLowFreqLists(BoxKey curkey, BoxDat& curdat);
    int SetupTreeHighFreqLists(BoxKey curkey, BoxDat& curdat);
    bool SetupTreeFind(BoxKey wntkey, BoxKey& reskey);
    bool SetupTreeAdjacent(BoxKey me, BoxKey yo);
    int RecursiveBoxInsert(std::queue< std::pair<BoxKey, BoxDat> >& tmpq,
                           bool first_pass);
    
    // AccLevel returns log_2(1 / epsilon), where epsilon is the target accuracy.
    int AccLevel();
    int SetupHighFreqCallLists();
    int GetExtPos();
    int GetHighFreqDirs();
    int SetupLowFreqOctree();
    int SetupLowFreqCallLists(); 

    // Travel up the octree and visit the boxes in the low frequency regime.
    // Compute outgoing non-directional equivalent densities using M2M.
    //
    // W is the width of the box
    // srcvec is the list of all boxes owned by this processor of width W
    // reqboxset is filled with all boxes whose information this processor needs
    // for the low-frequency downward pass
    int EvalUpwardLow(double W, std::vector<BoxKey>& srcvec,
                      std::set<BoxKey>& reqboxset);

    int EvalDownwardLow(double W, std::vector<BoxKey>& trgvec);
    int LowFreqUpwardPass(ldmap_t& ldmap, std::set<BoxKey>& reqboxset);
    int LowFreqDownwardComm(std::set<BoxKey>& reqboxset);
    int LowFreqDownwardPass(ldmap_t& ldmap);

    int HighFreqPass();

    int EvalUpwardHigh(double W, Index3 dir, std::vector<BoxKey>& srcvec);
    int EvalDownwardHigh(double W, Index3 dir, std::vector<BoxKey>& trgvec,
                         std::set<BoxAndDirKey>& affected_keys);

    int GatherLocalKeys();
    void ConstructLowFreqMap(ldmap_t& ldmap);
    int GatherDensities(ParVec<int,cpx,PtPrtn>& den);
    
    int UListCompute(BoxDat& trgdat);
    int VListCompute(double W, BoxDat& trgdat, Point3& trgctr, DblNumMat& uep,
		     DblNumMat& dcp, CpxNumVec& dnchkval,
                     NumTns<CpxNumTns>& ue2dc);
    int WListCompute(double W, BoxDat& trgdat, DblNumMat& uep);
    int XListCompute(BoxDat& trgdat, DblNumMat& dcp, DblNumMat& dnchkpos,
                     CpxNumVec& dnchkval);

    int LowFreqM2M(BoxKey& srckey, BoxDat& srcdat, DblNumMat& uep,
                   DblNumMat& ucp, NumVec<CpxNumMat>& uc2ue,
                   NumTns<CpxNumMat>& ue2uc);
    int LowFreqM2L(BoxKey& trgkey, BoxDat& trgdat, DblNumMat& dcp,
                   NumTns<CpxNumTns>& ue2dc, CpxNumVec& dneqnden, DblNumMat& uep,
                   NumVec<CpxNumMat>& dc2de);
    int LowFreqL2L(BoxKey& trgkey, BoxDat& trgdat, DblNumMat& dep,
                   NumTns<CpxNumMat>& de2dc, CpxNumVec& dneqnden);

    int HighFreqM2M(BoxAndDirKey& bndkey, NumVec<CpxNumMat>& uc2ue,
                    NumTns<CpxNumMat>& ue2uc);
    int HighFreqM2L(Index3 dir, BoxKey trgkey, DblNumMat& dcp, DblNumMat& uep);
    int HighFreqL2L(Index3 dir, BoxKey trgkey, NumVec<CpxNumMat>& dc2de,
		    NumTns<CpxNumMat>& de2dc,
                    std::set<BoxAndDirKey>& affected_keys);


    // Routines for communication
    int HighFreqInteractionListKeys(int level,
                                    std::set<BoxAndDirKey>& request_keys);

    int AllChildrenKeys(LevelBoxAndDirVec& vec,
                        std::vector<BoxAndDirKey>& req_keys);
    int HighFreqM2MLevelComm(int level);
    int HighFreqL2LLevelCommPre(int level);
    int HighFreqL2LLevelCommPost(int level, std::set<BoxAndDirKey>& affected_keys);
    int HighFreqM2LComm(int level,
                        std::set<BoxAndDirKey>& request_keys);
    int HighFreqL2LDataUp(BoxAndDirKey key, BoxAndDirDat& dat,
                          std::vector<int>& pids);
    static int HighFreqL2LDataUp_wrapper(BoxAndDirKey key, BoxAndDirDat& dat,
                                         std::vector<int>& pids);
    int HighFreqM2MDataUp(BoxAndDirKey key, BoxAndDirDat& dat,
                          std::vector<int>& pids);
    static int HighFreqM2MDataUp_wrapper(BoxAndDirKey key, BoxAndDirDat& dat,
                                         std::vector<int>& pids);

    // Tools for data distribution.
    void PrtnDirections(level_hdkeys_t& level_hdkeys,
                         std::vector<LevelBoxAndDirVec>& level_hf_vecs);
    int PrtnUnitLevel();
    int FormUnitPrtnMap(UnitLevelBoxPrtn& prtn, std::vector<int>& start_data,
                        std::vector<int>& end_data);

    int TransferBoxAndDirData(BoxAndDirKey key, BoxAndDirDat& dat,
                              std::vector<int>& pids);
    static int TransferBoxAndDirData_wrapper(BoxAndDirKey key, BoxAndDirDat& dat,
                                             std::vector<int>& pids);
    int TransferUnitLevelData(BoxKey key, BoxDat& dat,
                              std::vector<int>& pids);
    static int TransferUnitLevelData_wrapper(BoxKey key, BoxDat& dat,
                                             std::vector<int>& pids);
    int TransferDataToLevels();

    int CleanLevel(int level);
    int CleanBoxvec();
    void DeleteEmptyBoxes(std::map<BoxKey, BoxDat>& data);
};

//-------------------
int serialize(const PtPrtn&, std::ostream&, const std::vector<int>&);
int deserialize(PtPrtn&, std::istream&, const std::vector<int>&);
//-------------------
int serialize(const BoxKey&, std::ostream&, const std::vector<int>&);
int deserialize(BoxKey&, std::istream&, const std::vector<int>&);
//-------------------
int serialize(const BoxDat&, std::ostream&, const std::vector<int>&);
int deserialize(BoxDat&, std::istream&, const std::vector<int>&);
//-------------------
int serialize(const BoxAndDirDat&, std::ostream&, const std::vector<int>&);
int deserialize(BoxAndDirDat&, std::istream&, const std::vector<int>&);
//-------------------
int serialize(const BoxAndDirKey&, std::ostream&, const std::vector<int>&);
int deserialize(BoxAndDirKey&, std::istream&, const std::vector<int>&);

#endif  // _WAVE3D_HPP_
