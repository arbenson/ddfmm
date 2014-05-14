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
#include "DataCollection.hpp"
#include "file_io.h"
#include "numtns.hpp"
#include "parallel.hpp"
#include "wave3d.hpp"

#include <exception>
#include <iostream>
#include <sstream>
#include <string>

int optionsCreate(int argc, char** argv, std::map<std::string,
                  std::string>& options) {
    options.clear();
    for (int k = 1; k < argc; k += 2) {
      options[ std::string(argv[k]) ] = std::string(argv[k + 1]);
    }
    return 0;
}

std::string findOption(std::map<std::string, std::string>& opts,
                  std::string option) {
    std::map<std::string, std::string>::iterator mi = opts.find(option);
    if (mi == opts.end()) {
        std::cerr << "Missing option " << option << std::endl;
        return "";
    }
    return mi->second;
}

int main(int argc, char** argv) {
#ifndef RELEASE
    // When not in release mode, we catch all errors so that we can print the
    // manual call stack.
    try {
#endif
        MPI_Init(&argc, &argv);
        int mpirank, mpisize;
        getMPIInfo(&mpirank, &mpisize);

        Wave3d wave("wave3d_");

        //0. init and get options
        srand48(time(NULL));
        CHECK_TRUE(argc % 1 == 0);
        double t0, t1;
        std::vector<int> all(1,1);
        std::map<std::string, std::string> opts;
        optionsCreate(argc, argv, opts);
        std::map<std::string, std::string>::iterator mi;
        std::string opt;

        //1. read data
        opt = findOption(opts, "-posfile");
        if (opt.empty()) {
            return 0;
        }

        std::istringstream piss;
        SAFE_FUNC_EVAL( SeparateRead(opt, piss) );
        SAFE_FUNC_EVAL( deserialize(wave._positions, piss, all) );
        std::vector<int>& pos_breaks = wave._positions.prtn().ownerinfo();
	// numpts is the total number of discretization points.
        int numpts = pos_breaks.back();
	// If some processors were not assigned any data, we add them to the list.
	while (pos_breaks.size() < mpisize) {
	    pos_breaks.push_back(numpts + 1);
	}
        if (mpirank == 0) {
            std::cerr << "Total number of points: " << numpts << std::endl;
            std::cerr << "Done reading pos " << wave._positions.lclmap().size() << " "
                      << wave._positions.prtn().ownerinfo().size() << std::endl;
        }

        ParVec<int, cpx, PtPrtn> den;
        opt = findOption(opts, "-denfile");
        if (opt.empty()) {
            return 0;
        }
        std::istringstream diss;
        SAFE_FUNC_EVAL( SeparateRead(opt, diss) );
        SAFE_FUNC_EVAL( deserialize(den, diss, all) );
	// If some processors were not assigned any data, we add them to the list.
	std::vector<int>& den_breaks = den.prtn().ownerinfo();
	int last = den_breaks.back();
	while (den_breaks.size() < mpisize) {
	    den_breaks.push_back(last + 1);
	}

        if (mpirank == 0) {
            std::cerr << "Done reading den " << den.lclmap().size() << " "
                 << den.prtn().ownerinfo().size() << std::endl;
        }
        ParVec<int, cpx, PtPrtn> val;  // preset val to be the same as den

        val = den;
        if (mpirank == 0) {
            std::cerr << "Done setting val " << val.lclmap().size() << " "
                 << val.prtn().ownerinfo().size() << std::endl;
        }
        Kernel3d kernel(KERNEL_HELM);
        mi = opts.find("-kernel");
        // TODO (Austin): change this to the other format for getting an option
        if (mi != opts.end()) {
            std::istringstream ss(mi->second);
            int type;
            ss >> type;
            kernel.type() = type;
        }
        Mlib3d& mlib = wave._mlib;
        mlib.kernel() = kernel;
        mlib.NPQ() = 4;
        SAFE_FUNC_EVAL( mlib.setup(opts) );
        if (mpirank == 0) {
            std::cerr << "Done reading mlib " << mlib.w2ldmap().size() << " "
		      << mlib.w2hdmap().size() << std::endl;
        }

        double K;
        opt = findOption(opts, "-wave3d_K");
        if (opt.empty()) {
            return 0;
        }
        std::istringstream ss(opt);
        ss >> K;
        
        double NPW;
        opt = findOption(opts, "-wave3d_NPW");
        if (opt.empty()) {
            return 0;
        }
        std::istringstream ss2(opt);
        ss2 >> NPW;

        std::string geomfile;
        opt = findOption(opts, "-geomfile");
        if (opt.empty()) {
            return 0;
        }
        geomfile = opt;

        IntNumTns geomprtn;
        opt = findOption(opts, "-geomprtn");
        if (opt.empty()) {
          return 0;
        }
  
        std::istringstream giss;
        SAFE_FUNC_EVAL( SharedRead(opt, giss) );
        SAFE_FUNC_EVAL( deserialize(geomprtn, giss, all) );
        if (mpirank == 0) {
          std::cerr << "Done reading geomprtn " << geomprtn.m() << " "
                    << geomprtn.n() << " " << geomprtn.p() << std::endl
                    << geomprtn << std::endl;
        }

#if 0
        IntNumTns geomprtn;
        NewData(geomfile, K, NPW, mpisize, geomprtn);

        if (mpirank == 0) {
            std::cerr << "Done reading geomprtn "
                      << geomprtn.m() << " "
                      << geomprtn.n() << " " 
                      << geomprtn.p() << std::endl
                      << geomprtn << std::endl;
        }
#endif
        wave.kernel() = kernel;
        wave.geomprtn() = geomprtn;

        //2. setup
        t0 = MPI_Wtime();
        SAFE_FUNC_EVAL( wave.setup(opts) );
        t1 = MPI_Wtime();
        if (mpirank == 0) {
            std::cout << "wave setup used " << MPIDiffTime(t0, t1) << "secs " << std::endl;
        }

        // 3. eval
        double time_eval;
        t0 = MPI_Wtime();
        SAFE_FUNC_EVAL( wave.eval(den, val) );
        t1 = MPI_Wtime();
        if (mpirank == 0) {
            std::cout << "wave eval used " << MPIDiffTime(t0, t1) << "secs " << std::endl;
        }
        time_eval = MPIDiffTime(t0, t1);

        std::ostringstream oss;
        SAFE_FUNC_EVAL( serialize(val, oss, all) );

        opt = findOption(opts, "-valfile");
        if (opt.empty()) {
            return 0;
        }
        SAFE_FUNC_EVAL( SeparateWrite(opt, oss) );

        // 4. check
        IntNumVec chkkeyvec;
        opt = findOption(opts, "-chkfile");
        if (opt.empty()) {
            return 0;
        }
        std::istringstream iss;
        SAFE_FUNC_EVAL( SharedRead(opt, iss) );
        SAFE_FUNC_EVAL( deserialize(chkkeyvec, iss, all) );
        int numchk = chkkeyvec.m();
        t0 = MPI_Wtime();
        double relerr = wave.check(den, val, chkkeyvec);
        t1 = MPI_Wtime();
        if (mpirank == 0) {
            std::cout << "wave check used " << MPIDiffTime(t0, t1)
                      << "secs " << std::endl;
        }

        // 5. output results
        double time_drct = MPIDiffTime(t0, t1) * numpts / double(numchk);
        if (mpirank == 0) {
            printf("----------------------\n");
            printf("RESULT\n");
            printf("K  %.2e\n", wave.K());
            printf("Ta %.2e\n", time_eval);
            printf("Td %.2e\n", time_drct);
            printf("Rt %.2e\n", time_drct / time_eval);
            printf("Ea %.6e\n", relerr);
            printf("----------------------\n");
        }
        SAFE_FUNC_EVAL( MPI_Finalize() );
#ifndef RELEASE
    } catch( ... ) {
        int mpirank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
        std::cerr << "Process " << mpirank << " caught error." << std::endl;
        DumpCallStack();
    }
#endif
    return 0;
}
