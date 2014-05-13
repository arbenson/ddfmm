#ifndef _ACOU3D_HPP_
#define _ACOU3D_HPP_

#include "comobject.hpp"
#include "vec3t.hpp"
#include "numtns.hpp"
#include "kernel3d.hpp"
#include "wave3d.hpp"

#include <string>

using std::vector;
using std::pair;
using std::map;
using std::set;
using std::cerr;
using std::cout;
using std::ostream;
using std::istream;

class VertexData {
  Point3 pos;  // vertex position
  cpx den;     // density
};

//---------------------------------------------------------------------------
class Acou3d: public ComObject {
public:
  //input
  vector<Point3> _vertvec;
  vector<Index3> _facevec; // all triangle faces
  ParVec<int, VertexData, PtPrtn> _vertex_vec;
  Point3 _ctr;
  int _accu;
  Kernel3d _knlbie; // what bie kernel to use
  //local
  vector<double> _diavec; //diagonal (solid angle/4pi)
  vector<double> _arevec; //length of each segment
  vector<Point3> _posvec; //pos used in fmm
  vector<Point3> _norvec; //nor used in fmm
  Wave3d _wave;
  map<int,DblNumMat> _gauwgts;
  map<int,DblNumMat> _sigwgts;
  
  Acou3d(const std::string& p): ComObject(p), _wave(p+"_wave") {}
  ~Acou3d() {}
  int setup(vector<Point3>& vertvec, vector<Index3>& facevec,
	    Point3 ctr, int accu, Kernel3d knlbie);
  int eval(vector<Point3>& chk, vector<cpx>& den, vector<cpx>& val);
};

#endif
