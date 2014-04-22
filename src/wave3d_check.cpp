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
#include "wave3d.hpp"
#include "vecmatop.hpp"

//---------------------------------------------------------------------
int Wave3d::check(ParVec<int, cpx, PtPrtn>& den, ParVec<int, cpx, PtPrtn>& val,
                  IntNumVec& chkkeys, double& relerr) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::check");
#endif
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
  
    _self = this;
    int mpirank = getMPIRank();
    ParVec<int, Point3, PtPrtn>& pos = (*_posptr);
  
    // 1. get pos
    std::vector<int> all(1, 1);
    std::vector<int> chkkeyvec;
    for (int i = 0; i < chkkeys.m(); ++i) {
        chkkeyvec.push_back( chkkeys(i) );
    }
    pos.getBegin(chkkeyvec, all);
    pos.getEnd(all);
    std::vector<Point3> tmpsrcpos;
    for (std::map<int,Point3>::iterator mi = pos.lclmap().begin();
        mi != pos.lclmap().end(); ++mi) {
        if (pos.prtn().owner(mi->first) == mpirank) {
            tmpsrcpos.push_back(mi->second);
        }
    }
    std::vector<cpx> tmpsrcden;
    for (std::map<int,cpx>::iterator mi = den.lclmap().begin();
        mi != den.lclmap().end(); ++mi) {
        if (den.prtn().owner(mi->first) == mpirank) {
            tmpsrcden.push_back(mi->second);
        }
    }
  
    std::vector<Point3> tmptrgpos;
    for (int i = 0; i < chkkeyvec.size(); ++i) {
        tmptrgpos.push_back( pos.access(chkkeyvec[i]) );
    }

    DblNumMat srcpos(3, tmpsrcpos.size(), false, (double*)&(tmpsrcpos[0]));
    CpxNumVec srcden(tmpsrcden.size(), false, (cpx*)&(tmpsrcden[0]));
    DblNumMat trgpos(3, tmptrgpos.size(), false, (double*)&(tmptrgpos[0]));
    CpxNumVec trgval(tmptrgpos.size());

    CpxNumMat inter;
    SAFE_FUNC_EVAL( _kernel.kernel(trgpos, srcpos, srcpos, inter) );
    // If no points were assigned to this processor, then the trgval
    // should be zero.
    if (inter.n() != 0) {
	SAFE_FUNC_EVAL( zgemv(1.0, inter, srcden, 0.0, trgval) );
    } else {
	for (int i = 0; i < trgval.m(); ++i) {
	    trgval(i) = 0;
	}
    }

    CpxNumVec allval(trgval.m());
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    // Note: 2 doubles per complex number
    SAFE_FUNC_EVAL( MPI_Allreduce(trgval.data(), allval.data(), trgval.m() * 2,
				  MPI_DOUBLE,
                                  MPI_SUM, MPI_COMM_WORLD) );
  
    // 2. get val
    val.getBegin(chkkeyvec, all);
    val.getEnd(all);
    CpxNumVec truval(chkkeyvec.size());
    for (int i = 0; i < chkkeyvec.size(); ++i) {
        truval(i) = val.access(chkkeyvec[i]);
    }
  
    CpxNumVec errval(chkkeyvec.size());
    for (int i = 0; i < chkkeyvec.size(); ++i) {
        errval(i) = allval(i) - truval(i);
    }

    if (mpirank == 0) {
      std::cout << "computed: " << std::endl << truval 
		<< "actual: " << std::endl << allval << std::endl;

    }

    double tn = sqrt( energy(truval) );
    double en = sqrt( energy(errval) );
    relerr = en / tn;
  
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}
