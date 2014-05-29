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
    inter.resize(M, N);
    double TWO_PI = 2 * M_PI;
    cpx I(0, 1);
    double mindif2 = _mindif * _mindif;

    if (_type == KERNEL_HELM_SINGLE) {
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

        DblNumMat& kr = r;  // 2pi * r
        mat_dscale(M, N, kr, TWO_PI);
    
        DblNumMat skr(M, N), ckr(M, N);
        mat_dsincos(M, N, kr, skr, ckr);
    
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
    
        DblNumMat& kr = r;  // 2pi * r
        mat_dscale(M, N, kr, TWO_PI);
    
        DblNumMat skr(M, N), ckr(M, N);
        mat_dsincos(M, N, kr, skr, ckr);
    
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                inter(i, j) = cpx(ckr(i, j), skr(i, j));
            }
        }
    } else if (_type == KERNEL_HELM_MIXED) {
        // Mixed single + double layer potential.
        // D(x, y) = (exp^{i * 2pi * r} (1 - i * 2pi * r) / r^3) * (r \cdot n)
        // S(x, y) = exp^{i * 2pi * r} / r
        // Evaluate: D(x, y) - i * eta * S(x, y) 
        //           = (1 / r) * exp^{i * 2pi * r} * 
        //             ((r \cdot n)(1 - i * 2pi r) / r^2 - i * eta)
        DblNumMat rr(M, N);  // r^2
	DblNumMat rn(M, N);  // r \cdot n
	for (int j = 0; j < N; ++j) {
            for (int i = 0; i < M; ++i) {
                double x = trgpos(0, i) - srcpos(0, j);
                double y = trgpos(1, i) - srcpos(1, j);
                double z = trgpos(2, i) - srcpos(2, j);
		double n0 = srcnor(0, j);
		double n1 = srcnor(1, j);
		double n2 = srcnor(2, j);
		rn(i, j) = x * n0 + y * n1 + z * n2;
                rr(i, j) = x * x + y * y + z * z;
                if (rr(i, j) < mindif2) {
                    rr(i, j) = 1;
		    rn(i, j) = 0;
                }
	    }
        }

        DblNumMat r(M, N);
        mat_dsqrt(M, N, rr, r);
    
	DblNumMat ir2(M, N); // 1 / r^2
        mat_dinv(M, N, rr, ir2);

        DblNumMat ir(M, N);  // 1 / r
        mat_dinv(M, N, r, ir);

        DblNumMat kr = r;  // 2pi * r
        mat_dscale(M, N, kr, TWO_PI);

	// sin and cos of 2pi * r
        DblNumMat skr(M, N), ckr(M, N);
        mat_dsincos(M, N, kr, skr, ckr);

	double eta = TWO_PI;
	for (int j = 0; j < N; ++j) {
	    for (int i = 0; i < M; ++i) {
	        // exp^{i * 2pi * r} / r
                cpx exp = cpx(ckr(i, j), skr(i, j)) * ir(i, j);
	        // (r \cdot n) * (1 - i * 2pi * r) / r^2
	        cpx tmp1 = rn(i, j) * ir2(i, j) * (- cpx(1, kr(i, j)));
		// - i * eta
	        cpx tmp2 = cpx(0, eta);
	        inter(i, j) = exp * (tmp1 - tmp2);
            }
        }
    } else if (_type == KERNEL_HELM_DOUBLE) {
	for (int j = 0; j < N; ++j) {
	    for (int i = 0; i < M; ++i) {
		double r0 = trgpos(0, i) - srcpos(0, j);
		double r1 = trgpos(1, i) - srcpos(1, j);
		double r2 = trgpos(2, i) - srcpos(2, j);
		double n0 = srcnor(0, j);
		double n1 = srcnor(1, j);
		double n2 = srcnor(2, j);
		double rr = r0*r0 + r1*r1 + r2*r2;
		double rn = r0*n0 + r1*n1 + r2*n2;
		double r = sqrt(rr);
		double rrr = rr * r;
		if (r < _mindif) {
		    inter(i, j) = 0;
		} else {
		    double Kr = TWO_PI * r;
		    inter(i, j) = 1.0/(4*M_PI) * cpx(cos(Kr),sin(Kr)) * (1.0-I*Kr) / rrr * rn;
		}
	    }
	}
    } else {
	std::cerr << "Unknown kernel type " << _type << std::endl;
        throw new std::exception();
    }
    return 0;
}
