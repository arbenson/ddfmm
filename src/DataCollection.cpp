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

#include "commoninc.hpp"

#include "mpi.h"

#include <ctime>
#include <iostream>

double MPIDiffTime(double t0, double t1) { return (t1 - t0); }

ParData GatherParData(double t0, double t1) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::GatherParData");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    double diff = MPIDiffTime(t0, t1);
    double *rbuf = new double[mpisize];

    MPI_Gather((void *)&diff, 1, MPI_DOUBLE, rbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    ParData data = {0., 0., 0., 0.};
    if (mpirank == 0) {
	data.max = rbuf[0];
	data.min = rbuf[0];
        for (int i = 0; i < mpisize; i++) {
            data.mean += rbuf[i];
	    if (rbuf[i] > data.max) {
		data.max = rbuf[i];
	    }
	    if (rbuf[i] < data.min) {
		data.min = rbuf[i];
	    }
        }
        data.mean /= mpisize;

        for (int i = 0; i < mpisize; i++) {
            data.var += (rbuf[i] - data.mean) * (rbuf[i] - data.mean);
        }
        data.var /= (mpisize - 1);
    }

    delete[] rbuf;
    return data;
}

// TODO (Austin): Combine this with GatherParData
CommData GatherCommData(int amt) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::GatherCommData");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    int *rbuf = new int[mpisize];

    MPI_Gather((void *)&amt, 1, MPI_INT, rbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    CommData data = {0., 0., 0, 0, 0};
    if (mpirank == 0) {
	data.max = rbuf[0];
	data.min = rbuf[0];
        for (int i = 0; i < mpisize; i++) {
            data.total += rbuf[i];
	    if (rbuf[i] > data.max)
		data.max = rbuf[i];
	    if (rbuf[i] < data.min)
		data.min = rbuf[i];
        }
        data.mean = static_cast<double>(data.total) / mpisize;
        for (int i = 0; i < mpisize; i++) {
            double curr = static_cast<double>(rbuf[i]);
            data.var += (curr - data.mean) * (curr - data.mean);
        }
        data.var /= (mpisize - 1);
    }

    delete[] rbuf;
    return data;
}

void PrintParData(ParData data, std::string message) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::PrintParData");
#endif
    int mpirank = getMPIRank();
    if (mpirank == 0) {
	std::cout << message << std::endl
		  << "mean: " << data.mean << std::endl
		  << "var: "  << data.var  << std::endl
		  << "max: "  << data.max  << std::endl
		  << "min: "  << data.min  << std::endl;
    }
}

void PrintCommData(CommData data, std::string message) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::PrintCommData");
#endif
    int mpirank = getMPIRank();
    if (mpirank == 0) {
	std::cout << message << std::endl
	          << "total: " << data.total << std::endl
		  << "mean: " << data.mean << std::endl
		  << "var: "  << data.var  << std::endl
		  << "max: "  << data.max  << std::endl
		  << "min: "  << data.min  << std::endl;
    }
}

