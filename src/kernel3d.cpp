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
#include "kernel3d.hpp"
#include "vecmath.hpp"

double Kernel3d::_mindif = 1e-8;

//---------------------------------------------------------------------------
// exp(i*K*r)/r
// Computes 
int Kernel3d::kernel(const DblNumMat& trgpos, const DblNumMat& srcpos,
                     const DblNumMat& srcnor, CpxNumMat& inter) {
#ifndef RELEASE
    CallStackEntry entry("Kernel3d::kernel");
#endif
    int M = trgpos.n();
    int N = srcpos.n();
    double K = 2 * M_PI;
    cpx I(0, 1);
    double mindif2 = _mindif * _mindif;

    if (_type == KERNEL_HELM) {
        DblNumMat r2(M, N);
	for (int j = 0; j < N; ++j) {
            for (int i = 0; i < M; ++i) {
                double x = trgpos(0, i) - srcpos(0, j);
                double y = trgpos(1, i) - srcpos(1, j);
                double z = trgpos(2, i) - srcpos(2, j);
                double tmp = x * x + y * y + z * z;
                r2(i, j) = tmp;
                if (tmp < mindif2) {
                    r2(i, j) = 1;
                }
	  }
        }

        DblNumMat r(M, N);
        mat_dsqrt(M, N, r2, r);
    
        DblNumMat& ir = r2;  // 1 / r
        mat_dinv(M, N, r, ir);

        DblNumMat& kr = r;  // Kr
        mat_dscale(M, N, kr, K);
    
        DblNumMat skr(M, N), ckr(M, N);
        mat_dsincos(M, N, kr, skr, ckr);
    
        inter.resize(M, N);
	for (int j = 0; j < N; ++j) {
	    for (int i = 0; i < M; ++i) {
                inter(i, j) = cpx( ckr(i, j) * ir(i, j), skr(i, j) * ir(i, j) );
            }
        }
    } else if (_type == KERNEL_EXPR) {
        DblNumMat r2(M, N);
	for (int j = 0; j < N; ++j) {
	    for (int i = 0; i < M; ++i) {
                double x = trgpos(0, i) - srcpos(0, j);
                double y = trgpos(1, i) - srcpos(1, j);
                double z = trgpos(2, i) - srcpos(2, j);
                double tmp = x * x + y * y + z * z;
                r2(i, j) = tmp;
                if (tmp < mindif2) {
                    r2(i, j) = 1;
                }
            }
        }

        DblNumMat r(M, N);
        mat_dsqrt(M, N, r2, r);
    
        DblNumMat& kr = r;  // Kr
        mat_dscale(M, N, kr, K);
    
        DblNumMat skr(M, N), ckr(M, N);
        mat_dsincos(M, N, kr, skr, ckr);
    
        inter.resize(M, N);
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                inter(i, j) = cpx(ckr(i, j), skr(i, j));
            }
        }
    } else {
	std::cerr << "Unknown kernel type " << _type << std::endl;
        throw new std::exception();
    }
    return 0;
}
