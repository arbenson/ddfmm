#ifndef _KERNEL3D_HPP_
#define _KERNEL3D_HPP_

#include "nummat.hpp"
#include "vec3t.hpp"

using std::vector;

enum {
  KNL_HELM = 0,
  KNL_EXPR = 1
};


class Kernel3d
{
protected:
  static double _mindif;
  int _type;
public:
  Kernel3d(int t=KNL_HELM): _type(t) {;}
  ~Kernel3d() {;}
  int& type() { return _type; }
  
  int dim() const { return 3; }
  int sdof() const { return 1; }
  int tdof() const { return 1; }
  
  int kernel(const DblNumMat& trgpos, const DblNumMat& srcpos, const DblNumMat& srcnor, CpxNumMat& mat);
};

#endif
