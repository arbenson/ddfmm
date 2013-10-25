#include "kernel3d.hpp"
#include "vecmath.hpp"

double Kernel3d::_mindif = 1e-8;

//---------------------------------------------------------------------------
// exp(i*K*r)/r
int Kernel3d::kernel(const DblNumMat& trgpos, const DblNumMat& srcpos,
                     const DblNumMat& srcnor, CpxNumMat& inter)
{
#ifndef RELEASE
    CallStackEntry entry("Kernel3d::kernel");
#endif
    int M = trgpos.n();
    int N = srcpos.n();
    double K = 2*M_PI;
    cpx I(0,1);
    double mindif2 = _mindif * _mindif;

    if (_type == KNL_HELM) {
        //-------------------------------
        DblNumMat r2(M,N);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                double x = trgpos(0,i) - srcpos(0,j);
                double y = trgpos(1,i) - srcpos(1,j);
                double z = trgpos(2,i) - srcpos(2,j);
                double tmp = x*x + y*y + z*z;
                r2(i,j) = tmp;
                if (tmp < mindif2) {
                    r2(i,j) = 1;
                }
            }
        }

        DblNumMat r(M,N);
        mat_dsqrt(M, N, r2, r);
    
        DblNumMat& ir = r2;  //1/r
        mat_dinv(M, N, r, ir);

        DblNumMat& kr = r; //Kr
        mat_dscale(M, N, kr, K);
    
        DblNumMat skr(M,N), ckr(M,N);
        mat_dsincos(M, N, kr, skr, ckr);
    
        inter.resize(M,N);
        for (int i = 0; i < M; i++)     {
            for(int j=0; j<N; j++) {
                inter(i,j) = cpx( ckr(i,j)*ir(i,j), skr(i,j)*ir(i,j) );
            }
        }
    } else if (_type == KNL_EXPR) {
        //-------------------------------
        DblNumMat r2(M,N);
        for (int i=0; i<M; i++) {
            for (int j=0; j<N; j++) {
                double x = trgpos(0,i) - srcpos(0,j);
                double y = trgpos(1,i) - srcpos(1,j);
                double z = trgpos(2,i) - srcpos(2,j);
                double tmp = x*x + y*y + z*z;
                r2(i,j) = tmp;
                if (tmp < mindif2) {
                    r2(i,j) = 1;
                }
            }
        }

        DblNumMat r(M,N);
        mat_dsqrt(M, N, r2, r);
    
        DblNumMat& kr = r; //Kr
        mat_dscale(M, N, kr, K);
    
        DblNumMat skr(M,N), ckr(M,N);
        mat_dsincos(M, N, kr, skr, ckr);
    
        inter.resize(M,N);
        for(int i = 0; i < M; i++) {
            for(int j = 0; j < N; j++) {
                inter(i,j) = cpx( ckr(i,j), skr(i,j) );
            }
        }
    } else {
        //-------------------------------
	std::cerr << "Unknown kernel type " << _type << std::endl;
        iA(0);
    }
    return 0;
}
