/* Distributed Directional Fast Multipole Method
   Copyright (C) 2013 Austin Benson, Lexing Ying, and Jack Poulson

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
#ifndef _BLAS_H_
#define _BLAS_H_

#include <complex>
typedef std::complex<float> cpx8;
typedef std::complex<double> cpx16;

extern "C" {
    void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
		double *alpha, double *a, int *lda, double *b, int *ldb,
		double *beta, double *c, int *ldc);
    void dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda,
		double *x, int *incx, double *beta, double *y, int *incy);
    
    void zgemm_(char *transa, char *transb, int *m, int *n, int *k, cpx16 *alpha,
		cpx16 *a, int *lda, cpx16 *b, int *ldb, cpx16 *beta, cpx16 *c,
		int *ldc); 
    void zgemv_(char *trans, int *m, int *n, cpx16 *alpha, cpx16 *a, int *lda,
		cpx16 *x, int *incx, cpx16 *beta, cpx16 *y, int *incy);
    void dscal_(int* n, double* alpha, double* X, int* incr);
    void zdotc_(cpx16* result, int* n, cpx16 *x, int* incx, cpx16* y, int* incy);
}

#endif  // _BLAS_H_
