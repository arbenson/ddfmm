
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

// z = alpha x + beta y
void Add(CpxNumVec& x, CpxNumVec& y, CpxNumVec& z, cpx alpha, cpx beta) {
    CHECK_TRUE(x.m() == y.m() && x.m() == z.m() && x.m() > 0);
    for (int i = 0; i < x.m(); ++i) {
	z(i) = alpha * x(i) + beta * y(i);
    }
}

// Compute ||x||_2
double Norm2(CpxNumVec& x) {
    return sqrt(std::abs(InnerProduct(x, x)));
}

// Add in k-th column.
void UpdateQ(CpxNumMat& Q, CpxNumVec& r, cpx curr_resid, int k) {
    for (int i = 0; i < Q.m(); ++i) {
	Q(i, k) = r(i) / curr_resid;
    }
}

// Update H matrix in GMRES
void UpdateResidual(CpxNumMat& H, CpxNumMat& Q, CpxNumVec& r, int iter) {
    for (int i = 0; i < iter; ++i) {
	CpxNumVec qi(r.m(), false, Q.data() + i * r.m());
	H(i, iter) = InnerProduct(qi, r);
	Add(r, qi, r, cpx(1, 0), -H(i, iter));
    }
}
    
// Solve y_* = \arg\min_y \| h10 e_1 - H y \|_2
// x_* = x_0 + Q y
void FormSolution(CpxNumMat& H, CpxNumMat& Q, double h10,
		  CpxNumVec& x0, int num_iters) {
    // TODO: implement
}

// Very simple GMRES (following GVL).
void GMRES(CpxNumVec& b, CpxNumVec& x0,
	   std::function<void (CpxNumVec& x, CpxNumVec& y)> Apply,
	   double tol, int max_iter) {
    int mpirank = getMPIRank();
    // r = b - Ax
    CpxNumVec result(x0.m());
    CpxNumVec r(x0.m());
    CpxNumVec q(x0.m());
    Apply(x0, result);
    Add(b, result, r, cpx(1, 0), cpx(-1, 0));
    double h10 = Norm2(r);
    double normb = Norm2(b);

    CpxNumMat Q(x0.m(), max_iter);
    CpxNumMat H(max_iter + 1, max_iter);
    setvalue(H, cpx(0, 0));

    double curr_resid = h10;
    int iter = 0;
    while (curr_resid > tol * normb && iter < max_iter) {
	std::cout << "residual: " << curr_resid << std::endl;
	// q_{k+1} = r_k / h_{k+1, k}
	UpdateQ(Q, r, curr_resid, iter);

	// Update iteration
	++iter;

	// r_{k+1} = Aq_{k+1}
	Apply(q, r);

	// Get the new residual, fill in H
	UpdateResidual(H, Q, r, iter);
	
	// H_{k + 1, k} = || r_{k} ||
	curr_resid = Norm2(r);
	H(iter, iter - 1) = curr_resid;
    }

    // We are now satisfied with the residual.  Solve the problem.
    FormSolution(H, Q, h10, x0, iter);
}
