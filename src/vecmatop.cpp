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
#include "blas.h"
#include "numvec.hpp"
#include "nummat.hpp"
#include "numtns.hpp"
#include "vecmatop.hpp"

//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta,
	  DblNumMat& C) {
#ifndef RELEASE
    CallStackEntry entry("dgemm");
#endif
    CHECK_TRUE_MSG(A.m() == C.m(), "A.m != C.m");
    CHECK_TRUE_MSG(A.n() == B.m(), "A.n != B.m" );
    CHECK_TRUE_MSG(B.n() == C.n(), "B.n != C.n");
    SAFE_FUNC_EVAL( dgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
    return 0;
}
// ---------------------------------------------------------------------- 
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta,
	  double* C) {
#ifndef RELEASE
    CallStackEntry entry("dgemm");
#endif
    char transa = 'N';
    char transb = 'N';
    CHECK_TRUE_MSG(m != 0 && n != 0 && k != 0, "Zero dimensions");
    dgemm_(&transa, &transb, &m, &n, &k,
	   &alpha, A, &m, B, &k, &beta, C, &m);
    return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta,
	  DblNumVec& Y) {
#ifndef RELEASE
    CallStackEntry entry("dgemv");
#endif
    CHECK_TRUE_MSG(Y.m() == A.m(), "Y.m != A.m");
    CHECK_TRUE_MSG(A.n() == X.m(), "A.n != X.m");
    SAFE_FUNC_EVAL( dgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int dgemv(int m, int n, double alpha, double* A, double* X, double beta,
	  double* Y) {
#ifndef RELEASE
    CallStackEntry entry("dgemv");
#endif
  char trans = 'N';
  CHECK_TRUE_MSG(m != 0 && n != 0, "Zero dimensions");
  int incx = 1;
  int incy = 1;
  dgemv_(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
  return 0;
}




// ---------------------------------------------------------------------- 
int zgemm(cpx alpha, const CpxNumMat& A, const CpxNumMat& B, cpx beta,
	  CpxNumMat& C) {
#ifndef RELEASE
    CallStackEntry entry("zgemm");
#endif
    CHECK_TRUE_MSG(A.m() == C.m(), "A.m != C.m");
    CHECK_TRUE_MSG(A.n() == B.m(), "A.n != B.m");
    CHECK_TRUE_MSG(B.n() == C.n(), "B.n != C.n");
    SAFE_FUNC_EVAL( zgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  return 0;
}
// ---------------------------------------------------------------------- 
int zgemm(int m, int n, int k, cpx alpha, cpx* A, cpx* B, cpx beta, cpx* C) {
#ifndef RELEASE
    CallStackEntry entry("zgemm");
#endif
    char transa = 'N';
    char transb = 'N';
    CHECK_TRUE_MSG(m != 0 && n != 0 && k != 0, "Zero dimensions");
    zgemm_(&transa, &transb, &m, &n, &k,
	   &alpha, A, &m, B, &k, &beta, C, &m);
    return 0;
}
//Y <- alpha A X + beta Y
// ---------------------------------------------------------------------- 
int zgemv(cpx alpha, const CpxNumMat& A, const CpxNumVec& X, cpx beta,
	  CpxNumVec& Y) {
#ifndef RELEASE
    CallStackEntry entry("zgemv");
#endif
    CHECK_TRUE_MSG(Y.m() == A.m(), "Y.m != A.m");
    CHECK_TRUE_MSG(A.n() == X.m(), "A.n != X.m");
    SAFE_FUNC_EVAL( zgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
    return 0;
}
// ---------------------------------------------------------------------- 
int zgemv(int m, int n, cpx alpha, cpx* A, cpx* X, cpx beta, cpx* Y) {
#ifndef RELEASE
    CallStackEntry entry("zgemv");
#endif
    CHECK_TRUE_MSG(m > 0 && A != NULL && X != NULL, "Bad data or dimensions");
    char trans = 'N';
    CHECK_TRUE_MSG(m != 0 && n != 0, "Zero dimensions");
    int incx = 1;
    int incy = 1;
    zgemv_(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
    return 0;
}
