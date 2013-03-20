#include "kernel3d.hpp"

double Kernel3d::_mindif = 1e-8;
using std::polar;

/*
extern "C"
{
  void vdsqrt_(int* n, double*, double*);
  void vdinv_(int* n, double*, double*);
  void vdsincos_(int* n, double*, double*, double*);
}
*/

//---------------------------------------------------------------------------
//exp(i*K*r)/r
int Kernel3d::kernel(const DblNumMat& trgpos, const DblNumMat& srcpos, const DblNumMat& srcnor, CpxNumMat& inter)
{
  int M = trgpos.n();
  int N = srcpos.n();
  int sdof = this->sdof();
  int tdof = this->tdof();
  double K = 2*M_PI;
  cpx I(0,1);
  
  double mindif2 = _mindif*_mindif;
  int TTL = M*N;

  if(_type==KNL_HELM) {
    //-------------------------------
    DblNumMat r2(M,N);
    for(int i=0; i<M; i++) {
      for(int j=0; j<N; j++) {
	double x = trgpos(0,i) - srcpos(0,j);	  double y = trgpos(1,i) - srcpos(1,j);	  double z = trgpos(2,i) - srcpos(2,j);
	double tmp = x*x + y*y + z*z;
	r2(i,j) = tmp;
	if(tmp<mindif2)
	  r2(i,j) = 1;
      }
    }
    DblNumMat r(M,N);
    //vdsqrt_(&TTL, r2.data(), r.data());
    for(int i=0; i<M; i++)      for(int j=0; j<N; j++)	r(i,j) = sqrt(r2(i,j));
    
    DblNumMat& ir = r2;  //1/r
    //vdinv_( &TTL, r.data(), ir.data());
    for(int i=0; i<M; i++)      for(int j=0; j<N; j++)	ir(i,j) = 1.0/(r(i,j));

    DblNumMat& kr = r; //Kr
    for(int i=0; i<M; i++)	for(int j=0; j<N; j++)	  kr(i,j) *= K;
    
    DblNumMat skr(M,N), ckr(M,N);
    //vdsincos_(&TTL, kr.data(), skr.data(), ckr.data());
    for(int i=0; i<M; i++)      for(int j=0; j<N; j++)	sincos(kr(i,j), &skr(i,j), &ckr(i,j));
    
    inter.resize(M,N);
    for(int i=0; i<M; i++)	for(int j=0; j<N; j++)	  inter(i,j) = cpx( ckr(i,j)*ir(i,j), skr(i,j)*ir(i,j) );
  } else if(_type==KNL_EXPR) {
    //-------------------------------
    DblNumMat r2(M,N);
    for(int i=0; i<M; i++) {
      for(int j=0; j<N; j++) {
	double x = trgpos(0,i) - srcpos(0,j);	  double y = trgpos(1,i) - srcpos(1,j);	  double z = trgpos(2,i) - srcpos(2,j);
	double tmp = x*x + y*y + z*z;
	r2(i,j) = tmp;
	if(tmp<mindif2)
	  r2(i,j) = 1;
      }
    }
    DblNumMat r(M,N);
    //vdsqrt_(&TTL, r2.data(), r.data());
    for(int i=0; i<M; i++)      for(int j=0; j<N; j++)	r(i,j) = sqrt(r2(i,j));
    
    DblNumMat& kr = r; //Kr
    for(int i=0; i<M; i++)	for(int j=0; j<N; j++)	  kr(i,j) *= K;
    
    DblNumMat skr(M,N), ckr(M,N);
    //vdsincos_(&TTL, kr.data(), skr.data(), ckr.data());
    for(int i=0; i<M; i++)      for(int j=0; j<N; j++)	sincos(kr(i,j), &skr(i,j), &ckr(i,j));
    
    inter.resize(M,N);
    for(int i=0; i<M; i++)	for(int j=0; j<N; j++)	  inter(i,j) = cpx( ckr(i,j), skr(i,j) );
  } else {
    //-------------------------------
    iA(0);
  }
  return 0;
}

/*
//---------------------------------------------------------------------------
int Kernel3d::kernel(const DblNumMat& trgpos, const DblNumMat& srcpos, const DblNumMat& srcnor, CpxNumMat& inter)
{
  int sdof = this->sdof();
  int tdof = this->tdof();
  double K = 2*M_PI;
  cpx I(0,1);
  inter.resize(trgpos.n(), srcpos.n());  clear(inter);
  
  for(int i=0; i<trgpos.n(); i++) {
	for(int j=0; j<srcpos.n(); j++) {
	  double x = trgpos(0,i) - srcpos(0,j);	  double y = trgpos(1,i) - srcpos(1,j);	  double z = trgpos(2,i) - srcpos(2,j);
	  double r2 = x*x + y*y + z*z;
	  double r = sqrt(r2);
	  if(r<_mindif) {
		inter(i,j) = 0;
	  } else {	  //inter(i,j) = exp(IK*r) / r;
		//double Kr = K*r;
		//inter(i,j) = (cos(Kr) + I*sin(Kr))/r;		//inter(i,j) = 1/r;
		inter(i,j) = polar(1.0/r,K*r);
	  }
	}
  }
  return 0;
}
*/
