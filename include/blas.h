#ifndef _BLAS_H_
#define _BLAS_H_

#include <complex>
typedef std::complex<float> cpx8;
typedef std::complex<double> cpx16;

extern "C" {
  void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
              double *a, int *lda, double *b, int *ldb, double *beta, double *c,
              int *ldc);
  void dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda,
              double *x, int *incx, double *beta, double *y, int *incy);
  void zgemm_(char *transa, char *transb, int *m, int *n, int *k, cpx16 *alpha,
              cpx16 *a, int *lda, cpx16 *b, int *ldb, cpx16 *beta, cpx16 *c,
              int *ldc); 
  void zgemv_(char *trans, int *m, int *n, cpx16 *alpha, cpx16 *a, int *lda,
              cpx16 *x, int *incx, cpx16 *beta, cpx16 *y, int *incy); 
}

#endif  // _BLAS_H_


