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
void UpdateQ(CpxNumMat& Q, CpxNumVec& r, cpx curr_resid, int k) {
#ifndef RELEASE
    CallStackEntry entry("UpdateQ");
#endif
    for (int i = 0; i < Q.m(); ++i) {
	Q(i, k) = r(i) / curr_resid;
    }
}


// Update H matrix in GMRES
void UpdateResidual(CpxNumMat& H, CpxNumMat& Q, CpxNumVec& r, int iter) {
#ifndef RELEASE
    CallStackEntry entry("UpdateResidual");
#endif
    for (int i = 0; i < iter; ++i) {
	CpxNumVec qi(r.m(), false, Q.data() + i * r.m());
	H(i, iter) = InnerProduct(qi, r);
	Add(r, qi, r, cpx(1, 0), -H(i, iter));
    }
}

// After the GMRES iterations have completed, this forms the solution vector.
void FormSolution(CpxNumMat& H, CpxNumMat& Q, double h10,
		  CpxNumVec& x0, int num_iters) {
#ifndef RELEASE
    CallStackEntry entry("FormSolution");
#endif
    // Solve y_* = \arg\min_y \| h10 e_1 - H y \|_2
    //   zgels (TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)
    char trans = 'N';
    int M = num_iters;
    int N = num_iters - 1;
    int NRHS = 1;
    int lda = H.m();
    int ldb = M;
    int lwork = -1;
    cpx opt_work;
    CpxNumVec b(num_iters);
    setvalue(b, cpx(0, 0));
    b(0) = h10;
    int info;

    // Workspace query
    zgels_(&trans, &M, &N, &NRHS, H.data(), &lda, b.data(), &ldb, &opt_work, &lwork, &info);
    CHECK_TRUE_MSG(info == 0, "Bad least-squares call");
    // Actual least squares solve
    lwork = static_cast<int>(std::real(opt_work));
    cpx *work = new cpx[lwork];
    zgels_(&trans, &M, &N, &NRHS, H.data(), &lda, b.data(), &ldb, work, &lwork, &info);
    CHECK_TRUE_MSG(info == 0, "Bad least-squares call");

    // Now solve x_* = x_0 + Q y.  y is stored in the vector b.
    zgemv(M, N, 1, Q.data(), b.data(), 1, x0.data());
}


// Very simple GMRES (following GVL).
void GMRES(CpxNumVec& b, CpxNumVec& x0,
	   std::function<void (CpxNumVec& x, CpxNumVec& y)> Apply,
	   double tol, int max_iter) {
#ifndef RELEASE
    CallStackEntry entry("GMRES");
#endif
    int mpirank = getMPIRank();
    // r = b - Ax
    CpxNumVec result(x0.m());
    CpxNumVec r(x0.m());
    Apply(x0, result);
    Add(b, result, r, cpx(1, 0), cpx(-1, 0));
    double h10 = Norm2(r);
    double normb = Norm2(b);
    if (mpirank == 0) {
	std::cout << "norm of right-hand-side: " << normb << std::endl;
    }

    CpxNumMat Q(x0.m(), max_iter);
    CpxNumMat H(max_iter + 1, max_iter);
    setvalue(H, cpx(0, 0));

    double curr_resid = h10;
    int iter = 0;
    while (curr_resid > tol * normb && iter < max_iter) {
	if (mpirank == 0) {
	    std::cout << "residual: " << curr_resid << std::endl;
	}
	// q_{k+1} = r_k / h_{k+1, k}
	UpdateQ(Q, r, curr_resid, iter);

	// Update iteration
	++iter;

	// r_{k+1} = Aq_{k+1}
	CpxNumVec q(x0.m(), false, Q.data() + (iter - 1) * x0.m());
	Apply(q, r);

	// Get the new residual, fill in H
	UpdateResidual(H, Q, r, iter);
	
	// H_{k + 1, k} = || r_{k} ||
	curr_resid = Norm2(r);
	H(iter, iter - 1) = curr_resid;

	MPI_Barrier(MPI_COMM_WORLD);
    }

    if (mpirank == 0) {
	std::cout << "Final residual: " << curr_resid << std::endl;
	std::cout << "Number of iterations: " << iter << std::endl;
    }

    // We are now satisfied with the residual.  Solve the problem.
    FormSolution(H, Q, h10, x0, iter);
}
