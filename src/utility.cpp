#include "commoninc.hpp"

int getMPIRank() {
#ifndef RELEASE
    CallStackEntry entry("getMPIRank");
#endif
    int mpirank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    return mpirank;
}

int getMPISize() {
#ifndef RELEASE
    CallStackEntry entry("getMPISize");
#endif
    int mpisize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    return mpisize;
}

int getMPIInfo(int *mpirank, int *mpisize) {
#ifndef RELEASE
    CallStackEntry entry("getMPIInfo");
#endif
    *mpirank = getMPIRank();
    *mpisize = getMPISize();
    return 0;
}
