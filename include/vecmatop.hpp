#ifndef _VECMATOP_HPP_
#define _VECMATOP_HPP_

#include "nummat.hpp"

//--------------------------------------------------
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C);
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C);

int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y);
int dgemv(int m, int n, double alpha, double* A, double* X, double beta, double* Y);

//--------------------------------------------------
int zgemm(cpx alpha, const CpxNumMat& A, const CpxNumMat& B, cpx beta, CpxNumMat& C);
int zgemm(int m, int n, int k, cpx alpha, cpx* A, cpx* B, cpx beta, cpx* C);

int zgemv(cpx alpha, const CpxNumMat& A, const CpxNumVec& X, cpx beta, CpxNumVec& Y);
int zgemv(int m, int n, cpx alpha, cpx* A, cpx* X, cpx beta, cpx* Y);

#endif
