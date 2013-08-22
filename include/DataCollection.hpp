#ifndef _DATA_COLLECTION_HPP_
#define _DATA_COLLECTION_HPP_

#include "mpi.h"

typedef struct ParData {
    double mean;
    double var;
    double max;
    double min;
} ParData;

typedef struct CommData {
    double mean;
    double var;
    int max;
    int min;
    int total;
} CommData;


ParData GatherParData(time_t t0, time_t t1) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::GatherParData");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    double diff = difftime(t1, t0);
    double *rbuf = new double[mpisize];

    MPI_Gather((void *)&diff, 1, MPI_DOUBLE, rbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    ParData data = {0., 0., 0., 0.};
    if (mpirank == 0) {
	data.max = rbuf[0];
	data.min = rbuf[0];
        for (int i = 0; i < mpisize; i++) {
            data.mean += rbuf[i];
	    if (rbuf[i] > data.max)
		data.max = rbuf[i];
	    if (rbuf[i] < data.min)
		data.min = rbuf[i];
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
        cout << message << endl
	     << "mean: " << data.mean << endl
	     << "var: "  << data.var  << endl
	     << "max: "  << data.max  << endl
	     << "min: "  << data.min  << endl;
    }
}

void PrintCommData(CommData data, std::string message) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::PrintCommData");
#endif
    int mpirank = getMPIRank();
    if (mpirank == 0) {
        cout << message << endl
	     << "total: " << data.total << endl
	     << "mean: " << data.mean << endl
	     << "var: "  << data.var  << endl
	     << "max: "  << data.max  << endl
	     << "min: "  << data.min  << endl;
    }
}

#endif  // _DATA_COLLECTION_HPP_
