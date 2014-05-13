#include "acou3d.hpp"
#include "serialize.hpp"
#include "trmesh.hpp"

#include <fstream>

//-----------------------------------
int Acou3d::setup(vector<Point3>& vertvec, vector<Index3>& facevec,
		  Point3 ctr, int accu, Kernel3d knlbie) {
    _vertvec = vertvec;
    _facevec = facevec;
    _ctr = ctr;
    _accu = accu;
    _knlbie = knlbie;
    std::cerr << "type " << _knlbie.type() <<std::endl;
    // Compute the diagonal scaling.
    TrMesh trmesh("");
    std::cout << "setting up trmesh" << std::endl;
    SAFE_FUNC_EVAL( trmesh.setup(vertvec, facevec) );
    std::cout << "done setting up trmesh" << std::endl;
    std::cout << "computing interior" << std::endl;
    SAFE_FUNC_EVAL( trmesh.compute_interior(_diavec) );
    std::cout << "done computing interior" << std::endl;
    for(int k=0; k<_diavec.size(); k++) {
        _diavec[k] /= (4*M_PI);
    }
    std::cout << "computing area" << std::endl;
    SAFE_FUNC_EVAL( trmesh.compute_area(_arevec) );
    std::cout << "done computing area" << std::endl;
    // Load the quadrature weights
    vector<int> all(1,1);
    std::ifstream gin("gauwgts.bin");
    SAFE_FUNC_EVAL( deserialize(_gauwgts, gin, all) );
    std::cerr<<"gauwgts size "<<_gauwgts.size()<< std::endl;
    std::ifstream lin("sigwgts.bin");
    SAFE_FUNC_EVAL( deserialize(_sigwgts, lin, all) );
    std::cerr<<"sigwgts size "<<_sigwgts.size()<< std::endl;
    return 0;
}

//-----------------------------------
int Acou3d::eval(vector<Point3>& chk, vector<cpx>& den, vector<cpx>& val) {
  DblNumMat& gauwgt = _gauwgts[5];
  //
  int numgau = gauwgt.m();
  int NF = _facevec.size();
  _posvec.clear();
  _norvec.clear();
  vector<cpx> denvec;
  vector<cpx> valvec;

  for (int fi = 0; fi < _facevec.size(); ++fi) {
      Index3& face = _facevec[fi];
      // Get the three vertices of the face.
      Point3 pos0 = _vertvec[face(0)];
      Point3 pos1 = _vertvec[face(1)];
      Point3 pos2 = _vertvec[face(2)];

      double are = _arevec[fi];
      Point3 nor = cross(pos1 - pos0, pos2 - pos0);
      nor = nor / nor.l2();

      // Current density at the corners.
      cpx den0 = den[face(0)];
      cpx den1 = den[face(1)];
      cpx den2 = den[face(2)];

      for (int gi = 0; gi < numgau; ++gi) {
          double loc0 = gauwgt(gi, 0);
	  double loc1 = gauwgt(gi, 1);
	  double loc2 = gauwgt(gi, 2);
	  double wgt  = gauwgt(gi, 3);
	  _posvec.push_back(loc0 * pos0 + loc1 * pos1 + loc2 * pos2);
	  _norvec.push_back(nor);
	  denvec.push_back((loc0 * den0 + loc1 * den1 + loc2 * den2) * (are * wgt));
      }
  }
  for (Point3& point : chk) {
      _posvec.push_back(point);
      _norvec.push_back(point);
      denvec.push_back(cpx(0, 0));
  }
#if 0
  SAFE_FUNC_EVAL( _wave.setup(_posvec,_norvec, _ctr, _accu, _knlbie) );
#endif
  valvec.resize(denvec.size(), cpx(0, 0));
#if 0
  SAFE_FUNC_EVAL( _wave.eval(denvec, valvec) );
#endif
  //double relerr;    SAFE_FUNC_EVAL( _wave.check(denvec,valvec,4,relerr) );    std::cerr<<"relative error "<<relerr<< std::endl;
  val.resize(chk.size());
  for(int ci=0; ci<chk.size(); ci++) {
    val[ci] = valvec[ numgau*NF + ci ];
  }
  return 0;
}
