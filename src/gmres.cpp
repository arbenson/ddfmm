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

#include "vecmatop.hpp"

#include "commoninc.hpp"
#include "nummat.hpp"
#include "numvec.hpp"

#include <functional>

#include "blas.h"
#include "math.h"
#include "mpi.h"

// Compute x^H y
cpx InnerProduct(CpxNumVec& x, CpxNumVec& y) {
#ifndef RELEASE
    CallStackEntry entry("InnerProduct");
#endif
    CHECK_TRUE(x.m() == y.m() && x.m() > 0);
    cpx *data_x = x.data();
    cpx *data_y = y.data();
    int incx = 1;
    int incy = 1;
    int m = x.m();
    cpx result;
    zdotc_(&result, &m, data_x, &incx, data_y, &incy);

    cpx global_result;
    MPI_Allreduce(&result, &global_result, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return global_result;
}


// Compute ||x||_2
double Norm2(CpxNumVec& x) {
#ifndef RELEASE
    CallStackEntry entry("Norm2");
#endif
    return sqrt(std::abs(InnerProduct(x, x)));
}


// z = alpha x + beta y
void Add(CpxNumVec& x, CpxNumVec& y, CpxNumVec& z, cpx alpha, cpx beta) {
#ifndef RELEASE
    CallStackEntry entry("Add");
#endif
    CHECK_TRUE(x.m() == y.m() && x.m() == z.m() && x.m() > 0);
    for (int i = 0; i < x.m(); ++i) {
	z(i) = alpha * x(i) + beta * y(i);
    }
}


// Add in k-th column.
void UpdateV(CpxNumMat& V, CpxNumVec& r, double norm2_r, int k) {
#ifndef RELEASE
    CallStackEntry entry("UpdateV");
#endif
    CHECK_TRUE(V.m() == r.m());
    for (int i = 0; i < V.m(); ++i) {
	V(i, k) = r(i) / norm2_r;
    }
}


// Update H matrix in GMRES
void UpdateResidual(CpxNumMat& H, CpxNumMat& V, CpxNumVec& r, int k) {
#ifndef RELEASE
    CallStackEntry entry("UpdateResidual");
#endif
    for (int i = 0; i < k; ++i) {
	CpxNumVec qi(V.m(), false, V.data() + i * V.m());
	H(i, k - 1) = InnerProduct(qi, r);
	Add(r, qi, r, cpx(1, 0), -H(i, k - 1));
    }
}

// After the GMRES iterations have completed, this forms the solution vector.
void FormSolution(CpxNumMat& H, CpxNumMat& V, double h10,
		  CpxNumVec& x0, int iter) {
#ifndef RELEASE
    CallStackEntry entry("FormSolution");
#endif
    // Solve y_* = \arg\min_y \| h10 e_1 - H y \|_2
    char trans = 'N';
    int M = iter + 1;
    int N = iter;
    int NRHS = 1;
    int lda = H.m();
    CpxNumVec b(M);
    setvalue(b, cpx(0, 0));
    b(0) = h10;
    int ldb = M;

    int lwork = -1;
    cpx opt_work;
    int info;
    // Workspace query
    zgels_(&trans, &M, &N, &NRHS, H.data(), &lda, b.data(), &ldb,
	   &opt_work, &lwork, &info);
    CHECK_TRUE_MSG(info == 0, "Bad least-squares call");

    int mpirank = getMPIRank();
    if (mpirank == 0) {
	std::cout << "H: " << H << std::endl;
    }

    // Actual least squares solve
    lwork = static_cast<int>(std::real(opt_work));
    cpx *work = new cpx[lwork];
    zgels_(&trans, &M, &N, &NRHS, H.data(), &lda, b.data(), &ldb,
	   work, &lwork, &info);
    CHECK_TRUE_MSG(info == 0, "Bad least-squares call");
    delete [] work;

    // Now solve x_* = x_0 + V y.  y is stored in the vector b.
    int m = x0.m();
    int n = iter;
    lda = V.m();
    CHECK_TRUE(lda == m);
    cpx alpha(1, 0);
    int incx = 1;
    int incy = 1;
    cpx beta(1, 0);
    zgemv_(&trans, &m, &n, &alpha, V.data(), &lda, b.data(), &incx,
	   &beta, x0.data(), &incy);
}


// Very simple GMRES with restarts (following Saad).
void GMRES(CpxNumVec& b, CpxNumVec& x0,
	   std::function<void (CpxNumVec& x, CpxNumVec& y)> Apply,
	   double tol, int max_iter) {
#ifndef RELEASE
    CallStackEntry entry("GMRES");
#endif
    int mpirank = getMPIRank();

    double normb = Norm2(b);
    if (mpirank == 0) {
	std::cout << "norm of right-hand-side: " << normb << std::endl;
    }
    int total_iter = 0;
    int restart_size = 5;
    double beta;
    double break_tol = 1e-14;

    while (total_iter < max_iter) {
	// r = b - Ax
	CpxNumVec result(x0.m());
	CpxNumVec w(x0.m());
	Apply(x0, w);
	Add(b, w, w, cpx(1, 0), cpx(-1, 0));
	
	// beta = || r ||_2
	double beta = Norm2(w);
	if (mpirank == 0) { std::cout << "Starting residual: " << beta << std::endl; }
	if (beta < tol * normb) { break; }

	// Setup V with first column = r / beta
	CpxNumMat V(x0.m(), restart_size);
	setvalue(V, cpx(0, 0));
	UpdateV(V, w, beta, 0);

	// Initialize Hessenberg matrix
	CpxNumMat H(restart_size + 1, restart_size);
	setvalue(H, cpx(0, 0));

	int j;
	for (j = 0; j < restart_size; ++j) {
	    // w_{j} = Av_{j}
	    CpxNumVec v(V.m(), false, V.data() + j * V.m());
	    Apply(v, w);
	    
	    // Update Hessenberg matrix and w
	    for (int i = 0; i <= j; ++i) {
		CpxNumVec vi(V.m(), false, V.data() + i * V.m());
		H(i, j) = InnerProduct(w, vi);
		Add(w, vi, w, cpx(1, 0), -H(i, j));
	    }

	    H(j + 1, j) = Norm2(w);
	    if (std::abs(H(j + 1, j)) < break_tol) {
		if (mpirank == 0) {
		    std::cout << "Early stop" << std::endl;
		}
		break;
	    }

	    ++total_iter;

	    if (j < restart_size - 1) {
		// Generate next column
		UpdateV(V, w, std::abs(H(j + 1, j)), j + 1);
	    }

	    MPI_Barrier(MPI_COMM_WORLD);
	}
	    
	if (mpirank == 0) { std::cout << "Forming solution" << std::endl; }
	FormSolution(H, V, beta, x0, j);
    }

    if (mpirank == 0) {
	std::cout << "Number of iterations: " << total_iter << std::endl;
	std::cout << "soln: " << x0 << std::endl;
    }
}
