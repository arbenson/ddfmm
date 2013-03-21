#include "blas.h"
#include "lapack.h"
#include "numvec.hpp"
#include "nummat.hpp"
#include "numtns.hpp"

#include "vecmatop.hpp"

using std::cerr;

//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  iC( dgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C)
{
  char transa = 'N';
  char transb = 'N';
  assert(m!=0 && n!=0 && k!=0);
  dgemm_(&transa, &transb, &m, &n, &k,
		 &alpha, A, &m, B, &k, &beta, C, &m);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  iC( dgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dgemv(int m, int n, double alpha, double* A, double* X, double beta, double* Y)
{
  char trans = 'N';
  assert(m!=0 && n!=0);
  int incx = 1;
  int incy = 1;
  dgemv_(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
  return 0;
}




// ---------------------------------------------------------------------- 
int zgemm(cpx alpha, const CpxNumMat& A, const CpxNumMat& B, cpx beta, CpxNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  iC( zgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int zgemm(int m, int n, int k, cpx alpha, cpx* A, cpx* B, cpx beta, cpx* C)
{
  char transa = 'N';
  char transb = 'N';
  assert(m!=0 && n!=0 && k!=0);
  zgemm_(&transa, &transb, &m, &n, &k,
		 &alpha, A, &m, B, &k, &beta, C, &m);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int zgemv(cpx alpha, const CpxNumMat& A, const CpxNumVec& X, cpx beta, CpxNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  iC( zgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int zgemv(int m, int n, cpx alpha, cpx* A, cpx* X, cpx beta, cpx* Y)
{
  iA(m>0 && A!=NULL && X!=NULL);
  char trans = 'N';
  assert(m!=0 && n!=0);
  int incx = 1;
  int incy = 1;
  zgemv_(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
  return 0;
}
