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
#include "parUtils.h"

#include "mpi.h"

#include <vector>

namespace par {
template <>
class Mpi_datatype<HFBoxAndDirectionKey> {
public:
    static MPI_Datatype value() {
        static bool         first = true;
        static MPI_Datatype datatype;
        if (first) {
            first = false;
            MPI_Type_contiguous(sizeof(HFBoxAndDirectionKey), MPI_BYTE, &datatype);
            MPI_Type_commit(&datatype);
        }
        return datatype;
    }
};
}

std::vector<int> BoxAndDirection(HFBoxAndDirectionKey& key) {
    std::vector<int> data(6);
    // Direction information.
    data[0] = key.second[0];
    data[1] = key.second[1];
    data[2] = key.second[2];
    // Box index.
    data[3] = key.first.second[0];
    data[4] = key.first.second[1];
    data[5] = key.first.second[2];
}

void Wave3d::PartitionDirections(level_hdkeys_t& level_hdkeys_out,
                                 level_hdkeys_t& level_hdkeys_inc) {
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    // Figure out which level is the starting level.
    int local_start_level = 0;
    for (int i = 0; i < level_hdkeys_out.size(); ++i) {
        if (level_hdkeys_out[i].size() > 0) {
            local_start_level = i;
            break;
        }
    }

    // Make sure that my starting level agrees with everyone else.
    int global_start_level = 0;
    MPI_Allreduce(&local_start_level, &global_start_level, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    CHECK_TRUE(global_start_level == local_start_level);

    // Sort keys amongst processes.
    for (int i = global_start_level; i < unitlevel(); ++i) {
        if (mpirank == 0) {
            std::cerr << "Partitioning level: " << i << std::endl;
        }
        std::vector<HFBoxAndDirectionKey>& curr_level_keys = level_hdkeys_out[i];
        CHECK_TRUE(curr_level_keys.size() > 0);
        par::bitonicSort(curr_level_keys, MPI_COMM_WORLD);
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
        // Communicate starting keys for each processor.
        std::vector<int> first_data = BoxAndDirection(curr_level_keys[0]);
        std::vector<int> first_recv_buf(first_data.size() * mpisize);
        SAFE_FUNC_EVAL(MPI_Alltoall(&first_data[0], first_data.size(), MPI_INT,
                                    &first_recv_buf[0], first_recv_buf.size(),
                                    MPI_INT, MPI_COMM_WORLD));
#if 0
        // Communicate ending keys for each processor.
        std::vector<int> end_data = BoxAndDirection(curr_level_keys.back());
        std::vector<int> end_recv_buf(end_data.size() * mpisize);
        SAFE_FUNC_EVAL(MPI_Alltoall(&end_data[0], end_data.size(), MPI_INT,
                                    &end_recv_buf[0], end_recv_buf.size(),
                                    MPI_INT, MPI_COMM_WORLD));
#endif
    }

}
