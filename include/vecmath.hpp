#ifndef _VECMATH_
#define _VECMATH_
#include <math.h>

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
#ifdef MKL
  int TTL = M * N;
  vdsqrt_(&TTL, in.data(), out.data());
#else
  for(int i = 0; i < M; i++) {
    for(int j = 0; j < N; j++)	{
      out(i,j) = sqrt(in(i,j));
    }
  }
#endif
  return 0;
}

int mat_dinv(int M, int N, DblNumMat& in, DblNumMat& out)
{
#ifdef MKL
  int TTL = M * N;
  vdinv_(&TTL, in.data(), out.data());
  return;
#else
  for(int i = 0; i < M; i++) {
    for(int j = 0; j < N; j++)	{
      out(i,j) = 1.0 / in(i,j);
    }
  }
#endif
  return 0;
}

int mat_dsincos(int M, int N, DblNumMat& in, DblNumMat& out_sin,
                 DblNumMat& out_cos)
{
#ifdef MKL
  int TTL = M * N;
  vdsincos_(&TTL, in.data(), out_sin.data(), out_cos.data());
  return;
#else
  for(int i = 0; i < M; i++) {
    for(int j = 0; j < N; j++)	{
      sincos(in(i,j), &out_sin(i,j), &out_cos(i,j));
    }
  }
#endif
  return 0;
}

int mat_dscale(int M, int N, DblNumMat& in, double K)
{
  for(int i = 0; i < M; i++) {
    for(int j = 0; j < N; j++) {
      in(i,j) *= K;
    }
  }
  return 0;
}

#endif