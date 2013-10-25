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
#ifndef _VECMATH_
#define _VECMATH_
#include <math.h>
#include "blas.h"

/*
 * This file is a wrapper around vector operations.  If MKL is available, these functions
 * call the corresponding MKL functions.  Otherwise, the operations are done with
 * manual looping.
 */

#ifdef MKL
extern "C"
{
    void vdsqrt_(int* n, double*, double*);
    void vdinv_(int* n, double*, double*);
    void vdsincos_(int* n, double*, double*, double*);
}
#endif

int mat_dsqrt(int M, int N, DblNumMat& in, DblNumMat& out)
{
#ifndef RELEASE
    CallStackEntry entry("mat_dsqrt");
#endif
#ifdef MKL
    int TTL = M * N;
    vdsqrt_(&TTL, in.data(), out.data());
#else
    for(int i = 0; i < M; i++) {
        for(int j = 0; j < N; j++)      {
            out(i,j) = sqrt(in(i,j));
        }
    }
#endif
    return 0;
}

int mat_dinv(int M, int N, DblNumMat& in, DblNumMat& out)
{
#ifndef RELEASE
    CallStackEntry entry("mat_dinv");
#endif
#ifdef MKL
    int TTL = M * N;
    vdinv_(&TTL, in.data(), out.data());
#else
    for(int i = 0; i < M; i++) {
        for(int j = 0; j < N; j++)      {
            out(i,j) = 1.0 / in(i,j);
        }
    }
#endif
    return 0;
}

int mat_dsincos(int M, int N, DblNumMat& in, DblNumMat& out_sin,
                DblNumMat& out_cos)
{
#ifndef RELEASE
    CallStackEntry entry("mat_dsincos");
#endif
#ifdef MKL
    int TTL = M * N;
    vdsincos_(&TTL, in.data(), out_sin.data(), out_cos.data());
#else
    for(int i = 0; i < M; i++) {
        for(int j = 0; j < N; j++)      {
#ifdef OS_X
            out_sin(i, j) = sin(in(i, j));
            out_cos(i, j) = cos(in(i, j));
#else
            sincos(in(i,j), &out_sin(i,j), &out_cos(i,j));
#endif
        }
    }
#endif
    return 0;
}

int mat_dscale(int M, int N, DblNumMat& in, double K)
{
#ifndef RELEASE
    CallStackEntry entry("mat_dscale");
#endif
    int TTL = M * N;
    int incr = 1;
    dscal_(&TTL, &K, in.data(), &incr);
    return 0;
}

#endif
