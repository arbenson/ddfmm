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
class Mpi_datatype<BoxAndDirKey> {
public:
    static MPI_Datatype value() {
        static bool         first = true;
        static MPI_Datatype datatype;
        if (first) {
            first = false;
            MPI_Type_contiguous(sizeof(BoxAndDirKey), MPI_BYTE, &datatype);
            MPI_Type_commit(&datatype);
        }
        return datatype;
    }
};
}

void BoxAndDirection(BoxAndDirKey& key, std::vector<int>& out_key) {
#ifndef RELEASE
    CallStackEntry entry("BoxAndDirection");
#endif
    out_key.resize(6);

    // Direction information.
    out_key[0] = key._dir[0];
    out_key[1] = key._dir[1];
    out_key[2] = key._dir[2];
    // Box index.
    out_key[3] = key._boxkey.second[0];
    out_key[4] = key._boxkey.second[1];
    out_key[5] = key._boxkey.second[2];
}

void FormPartitionMap(HFBoxAndDirMap& map, std::vector<int>& start_data,
                      std::vector<int>& end_data, int level) {
    CHECK_TRUE(start_data.size() == end_data.size());
    CHECK_TRUE(start_data.size() % 6 == 0);
    std::vector<BoxAndDirKey>& part = map.partition_;
    // Keys are represented as 6 integers:
    //    (x, y, z) direction
    //    (x, y, z) box index
    for (int i = 0; i < start_data.size(); i += 6) {
        Index3 dir(start_data[i], start_data[i + 1], start_data[i + 2]);
        Index3 ind(start_data[i + 3], start_data[i + 4], start_data[i + 5]);
        BoxKey boxkey(level, ind);
        part.push_back(BoxAndDirKey(boxkey, dir));
    }
    
    // We only need the starting keys to determine the partition.  However,
    // we also store the ending keys for debugging.
    std::vector<BoxAndDirKey>& end_part = map.end_partition_;
    for (int i = 0; i < end_data.size(); i += 6) {
        Index3 dir(end_data[i], end_data[i + 1], end_data[i + 2]);
        Index3 ind(end_data[i + 3], end_data[i + 4], end_data[i + 5]);
        BoxKey boxkey(level, ind);
        end_part.push_back(BoxAndDirKey(boxkey, dir));
    }
}

void ScatterKeys(level_hdkeys_t& level_hdkeys) {
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    // Get the size of the keys on each 
    int my_size = level_hdkeys.size();
    std::vector<int> sizes(mpisize);
    SAFE_FUNC_EVAL( MPI_Alltoall(&my_size, 1, MPI_INT, &sizes[0], 1, MPI_INT,
                                 MPI_COMM_WORLD) );

    int my_index = -1;
    std::vector<int> empty_procs;
    std::vector<int> nonempty_procs;
    for (int p = 0; p < sizes.size(); ++p) {
       if (sizes[p] > 0) {
           nonempty_procs.push_back(p);
           if (p == mpirank) {
               my_index = nonempty_procs.size() - 1;
           }
       } else {
           empty_procs.push_back(p);
           if (p == mpirank) {
               my_index = empty_procs.size() - 1;
           }
       }
    }

    if (my_size == 0) {
        // Compute who is going to send me a key
        int my_sender = nonempty_procs[my_index % nonempty_procs.size()];
        int count = 0;
        for (int i = 0; i < empty_procs.size(); ++i) {
            if (nonempty_procs[i % nonempty_procs.size()] == my_sender) {
                ++count;
            }
        }
        int num_to_recv = sizes[my_sender] / (count + 1);
        // MPI_Irecv
    } else {
        // Compute to whom I am going to send data
        std::vector<int> dest_procs;
        for (int i = 0; i < empty_procs.size(); ++i) {
            if (nonempty_procs[i % nonempty_procs.size()] == mpirank) {
                dest_procs.push_back(nonempty_procs[i]);
            }
        }
        for (int i = 0; i < dest_procs.size(); ++i) {
            // MPI_Isend
        }
    }
}

void Wave3d::PartitionDirections(level_hdkeys_t& level_hdkeys_out,
                                 level_hdkeys_t& level_hdkeys_inc) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::PartitionDirections");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    ScatterKeys(level_hdkeys_out);
    ScatterKeys(level_hdkeys_inc);

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
        std::vector<BoxAndDirKey>& curr_level_keys = level_hdkeys_out[i];
        CHECK_TRUE(curr_level_keys.size() > 0);
        par::bitonicSort(curr_level_keys, MPI_COMM_WORLD);
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

        // Communicate starting keys for each processor.
        std::vector<int> start_data;
        BoxAndDirection(curr_level_keys[0], start_data);
        std::vector<int> first_recv_buf(start_data.size() * mpisize);
        SAFE_FUNC_EVAL(MPI_Alltoall((void *)&start_data[0], start_data.size(), MPI_INT,
                                    (void *)&first_recv_buf[0], start_data.size(),
                                    MPI_INT, MPI_COMM_WORLD));
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

        // Communicate ending keys for each processor.
        std::vector<int> end_data;
        BoxAndDirection(curr_level_keys.back(), end_data);
        std::vector<int> end_recv_buf(end_data.size() * mpisize);
        SAFE_FUNC_EVAL(MPI_Alltoall(&end_data[0], end_data.size(), MPI_INT,
                                    &end_recv_buf[0], end_data.size(),
                                    MPI_INT, MPI_COMM_WORLD));
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

        HFBoxAndDirMap map;
        FormPartitionMap(map, start_data, end_data, i);
    }
}
