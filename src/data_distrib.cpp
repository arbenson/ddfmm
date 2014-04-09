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

#include "external/bitonic.h"
#include "external/dtypes.h"

#include "mpi.h"

#include <vector>

#define BOX_AND_DIR_KEY_MPI_SIZE (6)

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
    out_key.resize(BOX_AND_DIR_KEY_MPI_SIZE);

    // Direction information.
    out_key[0] = key._dir[0];
    out_key[1] = key._dir[1];
    out_key[2] = key._dir[2];
    // Box index.
    out_key[3] = key._boxkey.second[0];
    out_key[4] = key._boxkey.second[1];
    out_key[5] = key._boxkey.second[2];
}

void FillKeyVector(std::vector<BoxAndDirKey>& keys, std::vector<int>& data,
                   int level) {
#ifndef RELEASE
    CallStackEntry entry("FillKeyVector");
#endif
   keys.clear();
   // Keys are represented as 6 integers:
   //    (x, y, z) direction
   //    (x, y, z) box index
   for (int i = 0; i < data.size(); i += BOX_AND_DIR_KEY_MPI_SIZE) {
        Index3 dir(data[i], data[i + 1], data[i + 2]);
        Index3 ind(data[i + 3], data[i + 4], data[i + 5]);
        BoxKey boxkey(level, ind);
        keys.push_back(BoxAndDirKey(boxkey, dir));
    }
}

void FormPartitionMap(BoxAndDirLevelPrtn& map, std::vector<int>& start_data,
                      std::vector<int>& end_data, int level) {
#ifndef RELEASE
    CallStackEntry entry("FormPartitionMap");
#endif
    CHECK_TRUE(start_data.size() == end_data.size());
    CHECK_TRUE(start_data.size() / getMPISize() == BOX_AND_DIR_KEY_MPI_SIZE);

    std::vector<BoxAndDirKey>& part = map.partition_;
    FillKeyVector(part, start_data, level);
    
    // We only need the starting keys to determine the partition.  However,
    // we also store the ending keys for debugging.
    std::vector<BoxAndDirKey>& end_part = map.end_partition_;
    FillKeyVector(end_part, end_data, level);
}

void ScatterKeys(std::vector<BoxAndDirKey>& keys, int level) {
#ifndef RELEASE
    CallStackEntry entry("ScatterKeys");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    // Get the size of the keys on each 
    int my_size = keys.size();
    std::vector<int> sizes(mpisize);
    SAFE_FUNC_EVAL( MPI_Allgather(&my_size, 1, MPI_INT, &sizes[0], 1, MPI_INT,
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

    if (mpirank == 1) {
      for (int i = 0; i < sizes.size(); ++i) {
        std::cout << sizes[i] << " ";
      }
      std::cout << std::endl;
    }

    // Distribute data more or less evenly.
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
        std::vector<int> buf(num_to_recv * BOX_AND_DIR_KEY_MPI_SIZE);
        MPI_Status status;
        MPI_Recv(&buf[0], buf.size(), MPI_INT, my_sender, 0, MPI_COMM_WORLD,
                 &status);
        FillKeyVector(keys, buf, level);
    } else {
        // Compute to whom I am going to send data
        std::vector<int> dest_procs;
        for (int i = 0; i < empty_procs.size(); ++i) {
            if (nonempty_procs[i % nonempty_procs.size()] == mpirank) {
                dest_procs.push_back(nonempty_procs[i]);
            }
        }
        int num_to_send = sizes[mpirank] / (dest_procs.size() + 1);
        if (dest_procs.size() > 0) {
            MPI_Request *reqs = new MPI_Request[dest_procs.size()];
            for (int i = 0; i < dest_procs.size(); ++i) {
              // Fill in data for this proc
              std::vector<int> buf(num_to_send * BOX_AND_DIR_KEY_MPI_SIZE);
              for (int j = 0; j < num_to_send; ++j) {
                  std::vector<int> curr_key;
                  BoxAndDirection(keys[i * num_to_send + j], curr_key);
                  CHECK_TRUE(curr_key.size() == BOX_AND_DIR_KEY_MPI_SIZE);
                  for (int k = 0; k < BOX_AND_DIR_KEY_MPI_SIZE; ++k) {
                      buf[j * BOX_AND_DIR_KEY_MPI_SIZE + k] = curr_key[k];
                  }
              }
              // TODO(arbenson): make this non-blocking
              MPI_Send(&buf[0], buf.size(), MPI_INT, dest_procs[i], 0,
                       MPI_COMM_WORLD);
            }
            // Remove keys from my list that are now on other processors.
            std::vector<BoxAndDirKey> keys_to_keep;
            for (int i = num_to_send * dest_procs.size(); i < keys.size(); ++i) {
                keys_to_keep.push_back(keys[i]);
            }
            keys.clear();
            keys.resize(keys_to_keep.size());
            for (int i = 0; i < keys_to_keep.size(); ++i) {
                keys[i] = keys_to_keep[i];
            }
        }
    }
}

void Wave3d::PartitionDirections(level_hdkeys_t& level_hdkeys_out,
                                 level_hdkeys_t& level_hdkeys_inc,
                                 std::vector<LevelBoxAndDirVec>& level_hf_vecs) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::PartitionDirections");
#endif
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
    for (int level = global_start_level; level < UnitLevel(); ++level) {
        if (mpirank == 0) {
            std::cerr << "Partitioning level: " << level << std::endl;
        }
        std::vector<BoxAndDirKey>& curr_level_keys = level_hdkeys_out[level];
        ScatterKeys(curr_level_keys, level);
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
        CHECK_TRUE(curr_level_keys.size() > 0);
        bitonicSort(curr_level_keys, MPI_COMM_WORLD);
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

        // Communicate starting keys for each processor.
        std::vector<int> start_data;
        BoxAndDirection(curr_level_keys[0], start_data);
        std::vector<int> start_recv_buf(start_data.size() * mpisize);
        SAFE_FUNC_EVAL(MPI_Allgather((void *)&start_data[0], start_data.size(), MPI_INT,
                                    (void *)&start_recv_buf[0], start_data.size(),
                                    MPI_INT, MPI_COMM_WORLD));
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

        // Communicate ending keys for each processor.
        std::vector<int> end_data;
        BoxAndDirection(curr_level_keys.back(), end_data);
        std::vector<int> end_recv_buf(end_data.size() * mpisize);
        SAFE_FUNC_EVAL(MPI_Allgather(&end_data[0], end_data.size(), MPI_INT,
                                     &end_recv_buf[0], end_data.size(),
                                     MPI_INT, MPI_COMM_WORLD));
        FormPartitionMap(level_hf_vecs[level].prtn(), start_recv_buf, end_recv_buf, level);
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

        // Build my ParVec for this level.
        for (int i = 0; i < curr_level_keys.size(); ++i) {
            BoxAndDirKey key = curr_level_keys[i];
            BoxAndDirDat dummy;
            level_hf_vecs[level].insert(key, dummy);
        }
    }

}
