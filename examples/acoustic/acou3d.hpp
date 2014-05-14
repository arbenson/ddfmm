#ifndef _ACOU3D_HPP_
#define _ACOU3D_HPP_

#include "comobject.hpp"
#include "vec3t.hpp"
#include "numtns.hpp"
#include "kernel3d.hpp"
#include "wave3d.hpp"

#include <string>

class VertexData {
  Point3 pos;  // vertex position
  cpx den;     // density
};

//---------------------------------------------------------------------------
class Acoustic3D: public ComObject {
public:
  //input
  std::vector<Point3> _vertvec;
  std::vector<Index3> _facevec; // all triangle faces
  Point3 _ctr;
  int _accu;
  Kernel3d _knlbie; // what bie kernel to use
  //local
  std::vector<double> _diavec; //diagonal (solid angle/4pi)
  std::vector<double> _arevec; //length of each segment
  std::vector<Point3> _posvec; //pos used in fmm
  std::vector<Point3> _norvec; //nor used in fmm
  Wave3d _wave;
  map<int, DblNumMat> _gauwgts;
  map<int, DblNumMat> _sigwgts;
  
  Acoustic3D(const std::string& p): ComObject(p), _wave(p+"_wave") {}
  ~Acoustic3D() {}
  int setup(std::vector<Point3>& vertvec, std::vector<Index3>& facevec,
	    Point3 ctr, int accu, Kernel3d knlbie);
  int eval(std::vector<Point3>& chk, std::vector<cpx>& den, std::vector<cpx>& val);
};

#endif
