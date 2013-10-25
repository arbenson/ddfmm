/* Distributed Directional Fast Multipole Method
   Copyright (C) 2013 Austin Benson, Lexing Ying, and Jack Poulson

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
