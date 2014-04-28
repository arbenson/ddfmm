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
#include "parvec.hpp"

#include "external/sorting/bitonic.h"
#include "external/sorting/dtypes.h"

#include "mpi.h"

#include <vector>

#define BOX_KEY_MPI_SIZE (3)

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

template <>
class Mpi_datatype<BoxKey> {
public:
    static MPI_Datatype value() {
        static bool         first = true;
        static MPI_Datatype datatype;
        if (first) {
            first = false;
            MPI_Type_contiguous(sizeof(BoxKey), MPI_BYTE, &datatype);
            MPI_Type_commit(&datatype);
        }
        return datatype;
    }
};
}


void FormPrtnMap(BoxAndDirLevelPrtn& map, std::vector<BoxAndDirKey>& start_data,
                 std::vector<BoxAndDirKey>& end_data, int level) {
#ifndef RELEASE
    CallStackEntry entry("FormPrtnMap");
#endif
    CHECK_TRUE(start_data.size() == end_data.size());
    map.partition_ = start_data;
    // We only need the starting keys to determine the partition.  However,
    // we also store the ending keys for debugging.
    map.end_partition_ = end_data;
}

void AdvancePosInd(std::vector<int>& diffs, int& pos) {
    while (pos < static_cast<int>(diffs.size()) && diffs[pos] <= 0) { ++pos; }
}

void AdvanceNegInd(std::vector<int>& diffs, int& neg) {
    while (neg < static_cast<int>(diffs.size()) && diffs[neg] >= 0) { ++neg; }
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
    if (mpirank == 0) {
        std::cout << "sizes: " << std::endl;
        for (int i = 0; i < static_cast<int>(sizes.size()); ++i) {
            std::cout << i << ": " << sizes[i] << std::endl;
        }
    }

    unsigned int total_size = 0;
    for (int i = 0; i < static_cast<int>(sizes.size()); ++i) {
        total_size += static_cast<unsigned int>(sizes[i]);
    }
    unsigned int avg_size = total_size / mpisize;
    if (mpirank == 0) {
        std::cout << "Targeting " << avg_size << " keys per process." << std::endl;
    }
    CHECK_TRUE_MSG(avg_size > 0, "Not enough keys!");
    std::vector<int> diffs(mpisize);
    for (int i = 0; i < static_cast<int>(diffs.size()); ++i) {
        diffs[i] = sizes[i] - avg_size;
    }

    // Get starting indices of negative and positive positions
    int pos = 0;
    AdvancePosInd(diffs, pos);
    int neg = 0;
    AdvanceNegInd(diffs, neg);

    // amt_to_send[rank] is a list of (proc, amt) data to send
    std::vector< std::vector< std::pair<int, int> > > amt_to_send(mpisize);
    while (neg < diffs.size()) {
        if (-diffs[neg] <= diffs[pos]) {
            // positive can fill negative
            amt_to_send[pos].push_back(std::pair<int, int>(neg, -diffs[neg]));
            diffs[pos] += diffs[neg];  // Should _decrease_ diffs[neg]
            diffs[neg] = 0;
            AdvanceNegInd(diffs, neg);
            AdvancePosInd(diffs, pos);
        } else {
            // positive cannot fully fill negative
            amt_to_send[pos].push_back(std::pair<int, int>(neg, diffs[pos]));
            diffs[neg] += diffs[pos];
            diffs[pos] = 0;
            AdvancePosInd(diffs, pos);
        }
        CHECK_TRUE_MSG(!(neg < diffs.size() && pos >= diffs.size()),
                       "Ran out of keys to distribute.");
    }
    
    std::vector<MPI_Request> reqs;
    // Make all receive requests on this process.
    std::vector< std::vector<BoxAndDirKey> > recv_bufs(mpisize);
    for (int i = 0; i < mpisize; ++i) {
        for (int j = 0; j < static_cast<int>(amt_to_send[i].size()); ++j) {
            if (amt_to_send[i][j].first == mpirank) {
                recv_bufs[i].resize(amt_to_send[i][j].second);
                reqs.emplace_back();
                CHECK_TRUE_MSG(i != mpirank, "Trying to receive data to myself.");
                MPI_Irecv(&(recv_bufs[i][0]), recv_bufs[i].size(),
                          par::Mpi_datatype<BoxAndDirKey>::value(), i,
                          0, MPI_COMM_WORLD, &reqs.back());
            }
        }
    }

    // Make all send requests on this process.
    int key_ind = 0;
    std::vector< std::pair<int, int> >& my_sends = amt_to_send[mpirank];
    for (int i = 0; i < static_cast<int>(my_sends.size()); ++i) {
        reqs.emplace_back();
        CHECK_TRUE_MSG(my_sends[i].first != mpirank,
                       "Trying to send data to myself.");
        MPI_Isend(&keys[key_ind], my_sends[i].second, 
                  par::Mpi_datatype<BoxAndDirKey>::value(), my_sends[i].first,
                  0, MPI_COMM_WORLD, &reqs.back());
        key_ind += my_sends[i].second;
    }

    // Get the tail end of the keys that didn't get transferred.
    std::vector<BoxAndDirKey> keys_to_keep;
    for (int i = key_ind; i < static_cast<int>(keys.size()); ++i) {
        keys_to_keep.push_back(keys[i]);
    }

    // Wait for all of the data to transfer
    std::vector<MPI_Status> stats(reqs.size());
    MPI_Waitall(reqs.size(), &reqs[0], &stats[0]);
    keys.clear();  // all data has been handled 
    
    // Insert back into keys
    for (int i = 0; i < mpisize; ++i) {
        for (int j = 0; j < static_cast<int>(recv_bufs[i].size()); ++j) {
            keys.push_back(recv_bufs[i][j]);
        }
        std::vector<BoxAndDirKey>().swap(recv_bufs[i]);
    }
    for (int i = 0; i < static_cast<int>(keys_to_keep.size()); ++i) {
        keys.push_back(keys_to_keep[i]);
    }
    keys_to_keep.clear();
}

void Wave3d::PrtnDirections(level_hdkeys_t& level_hdkeys,
                            std::vector<LevelBoxAndDirVec>& level_hf_vecs) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::PrtnDirections");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    // Figure out which level is the starting level.
    int local_start_level = level_hdkeys.size();
    for (int i = 0; i < static_cast<int>(level_hdkeys.size()); ++i) {
        if (level_hdkeys[i].size() > 0) {
            local_start_level = i;
            break;
        }
    }

    // Get the starting level.
    int global_start_level = 0;
    MPI_Allreduce(&local_start_level, &global_start_level, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    // If we haven't determined the starting level yet, set it.
    if (_starting_level == 0) {
        _starting_level = global_start_level;
    }

    // Sort keys amongst processes.
    for (int level = global_start_level; level < UnitLevel(); ++level) {
        if (mpirank == 0) {
            std::cerr << "Partitioning level: " << level << std::endl;
        }
        std::vector<BoxAndDirKey>& curr_level_keys = level_hdkeys[level];
        ScatterKeys(curr_level_keys, level);
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
        CHECK_TRUE_MSG(curr_level_keys.size() > 0, "Empty keys");
        bitonicSort(curr_level_keys, MPI_COMM_WORLD);
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

        // Communicate starting keys for each processor.
        std::vector<BoxAndDirKey> start_recv_buf(mpisize);
        SAFE_FUNC_EVAL(MPI_Allgather(&curr_level_keys[0], 1,
                                     par::Mpi_datatype<BoxAndDirKey>::value(),
                                     &start_recv_buf[0], 1,
                                     par::Mpi_datatype<BoxAndDirKey>::value(),
                                     MPI_COMM_WORLD));
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

        // Communicate ending keys for each processor.
        std::vector<BoxAndDirKey> end_recv_buf(mpisize);
        SAFE_FUNC_EVAL(MPI_Allgather(&curr_level_keys[curr_level_keys.size() - 1], 1,
                                     par::Mpi_datatype<BoxAndDirKey>::value(),
                                     &end_recv_buf[0], 1,
                                     par::Mpi_datatype<BoxAndDirKey>::value(),
                                     MPI_COMM_WORLD));
        FormPrtnMap(level_hf_vecs[level].prtn(), start_recv_buf, end_recv_buf, level);
        SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

        // Build my ParVec for this level.
        for (int i = 0; i < static_cast<int>(curr_level_keys.size()); ++i) {
            BoxAndDirKey key = curr_level_keys[i];
            BoxAndDirDat dummy;
            level_hf_vecs[level].insert(key, dummy);
        }
    }

}

int Wave3d::FormUnitPrtnMap(UnitLevelBoxPrtn& prtn,
                            std::vector<int>& start_data,
                            std::vector<int>& end_data) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::FormUnitPrtnMap");
#endif
    int mpisize = getMPISize();
    CHECK_TRUE(start_data.size() == end_data.size());
    CHECK_TRUE(static_cast<int>(start_data.size()) / BOX_KEY_MPI_SIZE == mpisize);
    std::vector<BoxKey>& part = prtn.partition_;
    part.resize(0);
    std::vector<BoxKey>& end_part = prtn.end_partition_;
    end_part.resize(0);
    for (int i = 0; i < static_cast<int>(start_data.size());
       i += BOX_KEY_MPI_SIZE) {
        Index3 start_ind(start_data[i], start_data[i + 1], start_data[i + 2]);
        Index3 end_ind(end_data[i], end_data[i + 1], end_data[i + 2]);
        part.push_back(BoxKey(UnitLevel(), start_ind));
        end_part.push_back(BoxKey(UnitLevel(), end_ind));
    }
}


int Wave3d::PrtnUnitLevel() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::PrtnUnitLevel");
#endif
    std::vector<BoxAndDirKey>& keys_out = _level_prtns._hdkeys_out[UnitLevel()];
    std::vector<BoxAndDirKey>& keys_inc = _level_prtns._hdkeys_inc[UnitLevel()];

    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    // Start by distributing the keys more or less uniformly.
    ScatterKeys(keys_out, UnitLevel());
    ScatterKeys(keys_inc, UnitLevel());

    // TODO(arbenson): Use morton ordering here.
    // Deal with just the set of boxes.
    std::set<BoxKey> boxes_set;
    for (int i = 0; i < static_cast<int>(keys_out.size()); ++i) {
        boxes_set.insert(keys_out[i]._boxkey);
    }
    for (int i = 0; i < static_cast<int>(keys_inc.size()); ++i) {
        boxes_set.insert(keys_inc[i]._boxkey);
    }
    std::vector<BoxKey> boxes;
    boxes.insert(boxes.begin(), boxes_set.begin(), boxes_set.end());
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

    // Sort
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    CHECK_TRUE_MSG(boxes.size() > 0, "empty boxes before sort");
    bitonicSort(boxes, MPI_COMM_WORLD);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

    // Processor i sends starting box to processor i - 1.
    // Processor i receives starting box to processor i + 1.
    CHECK_TRUE_MSG(boxes.size() > 0, "empty boxes after sort");
    Index3 start_box = boxes[0]._index;
    Index3 nbr_start_box(0, 0, 0);
    MPI_Status status;
    MPI_Sendrecv(start_box.array(), BOX_KEY_MPI_SIZE, MPI_INT,
                 mpirank == 0 ? mpisize - 1 : mpirank - 1, 0,
                 nbr_start_box.array(), BOX_KEY_MPI_SIZE, MPI_INT, (mpirank + 1) % mpisize,
                 0, MPI_COMM_WORLD, &status);

    // Local adjustments for ending position if they happen to overlap.
    if (mpirank > 0) {
        while (boxes.size() > 0 && nbr_start_box == boxes.back()._index) {
            boxes.pop_back();
        }
    }
    CHECK_TRUE_MSG(boxes.size() > 0, "empty boxes after adjustments");

    // Communicate starting keys for each processor.
    std::vector<int> start_recv_buf(BOX_KEY_MPI_SIZE * mpisize);
    SAFE_FUNC_EVAL(MPI_Allgather(start_box.array(), BOX_KEY_MPI_SIZE, MPI_INT,
                                 &start_recv_buf[0], BOX_KEY_MPI_SIZE,
                                 MPI_INT, MPI_COMM_WORLD));
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
  
    // Communicate ending keys for each processor.
    std::vector<int> end_recv_buf(BOX_KEY_MPI_SIZE * mpisize);
    SAFE_FUNC_EVAL(MPI_Allgather(boxes.back()._index.array(), BOX_KEY_MPI_SIZE, MPI_INT,
                                 &end_recv_buf[0], BOX_KEY_MPI_SIZE,
                                 MPI_INT, MPI_COMM_WORLD));

    // Form the parvec
    FormUnitPrtnMap(_level_prtns._unit_vec.prtn(), start_recv_buf, end_recv_buf);
    // Copy to low-frequency boxvec partition.
    
    LowFreqBoxPrtn& prtn = _level_prtns._lf_boxvec.prtn();
    UnitLevelBoxPrtn& unit_prtn = _level_prtns._unit_vec.prtn();
    prtn.partition_ = unit_prtn.partition_;
    prtn.end_partition_ = unit_prtn.end_partition_;
    prtn.unit_level_ = UnitLevel();
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );

    // Note: parvec gets filled in later with call to TransferDataToLevels()
    return 0;  
}


// This transfer holds two purposes:
//    1. Pass the high-frequency interaction lists to the appropriate processors
//    2. Get the unit-level directions (we currently only have the boxes)
int Wave3d::TransferBoxAndDirData(BoxAndDirKey key, BoxAndDirDat& dat,
                                  std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::TransferBoxAndDirData");
#endif
    pids.clear();
    int owner = _level_prtns.Owner(key, false);
    int level = key._boxkey._level;
    if (level == UnitLevel()) {
        pids.push_back(owner);
    } else if (dat.interactionlist().size() > 0) {
        // It is an incoming direction iff the interaction list is nonempty
        pids.push_back(owner);
    }
    return 0;
}

int Wave3d::TransferBoxAndDirData_wrapper(BoxAndDirKey key, BoxAndDirDat& dat,
                                          std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::TransferBoxAndDirData_wrapper");
#endif
    return (Wave3d::_self)->TransferBoxAndDirData(key, dat, pids);
}

int Wave3d::TransferDataToLevels() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::TransferDataToLevels");
#endif
    int mpirank = getMPIRank();
    for (std::map<BoxAndDirKey, BoxAndDirDat>::iterator mi = _bndvec.lclmap().begin();
       mi != _bndvec.lclmap().end(); ++mi) {
        BoxAndDirKey key = mi->first;
        int owner = _level_prtns.Owner(key, false);
        if (owner != mpirank) {
            continue;
        }
        BoxAndDirDat dat = mi->second;
        int level = key._boxkey._level;
        if (level == UnitLevel()) {
            // Put unit-level directions into the appropriate vector
            _level_prtns._unit_vec.lclmap()[key] = dat;

        } else if (dat.interactionlist().size() > 0) {
            // Put high-frequency directions into the appropriate vector.
            // It is an incoming direction iff the interaction list is nonempty.
            LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_inc[level];
            vec.lclmap()[key] = dat;
        }
    }
    return 0;
}


// This transfer moves the unit level box data to the new partition.
int Wave3d::TransferUnitLevelData(BoxKey key, BoxDat& dat,
                                  std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::TransferUnitLevelData");
#endif
    pids.clear();
    if (key._level == UnitLevel()) {
        Index3 dummy_dir(1, 1, 1);
        BoxAndDirKey new_key(key, dummy_dir);
        pids.push_back(_level_prtns.Owner(new_key, false));
    }
    return 0;
}

int Wave3d::TransferUnitLevelData_wrapper(BoxKey key, BoxDat& dat,
                                          std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::TransferUnitLevelData_wrapper");
#endif
    return (Wave3d::_self)->TransferUnitLevelData(key, dat, pids);
}

void LevelPartitions::Init(int K) {
#ifndef RELEASE
    CallStackEntry entry("LevelPartitions::init");
#endif
    unit_level_ = static_cast<int>(round(log(K) / log(2)));
    _hf_vecs_out.resize(unit_level_ + 1);
    _hf_vecs_inc.resize(unit_level_ + 1);
    _hdkeys_out.resize(unit_level_ + 1);
    _hdkeys_inc.resize(unit_level_ + 1);
    _level_hdmap_out.resize(unit_level_ + 1);
    _level_hdmap_inc.resize(unit_level_ + 1);
}

void InsertIntoDirMap(std::map<Index3, std::vector<BoxKey> >& dir_map,
                      std::map<BoxAndDirKey, BoxAndDirDat>& curr_map) {
#ifndef RELEASE
    CallStackEntry entry("InsertIntoDirMap");
#endif
    for (std::map<BoxAndDirKey, BoxAndDirDat>::iterator mi = curr_map.begin();
         mi != curr_map.end(); ++mi) {
        BoxAndDirKey key = mi->first;
        Index3 dir = key._dir;
        BoxKey boxkey = key._boxkey;
        // TODO(arbenson): check to make sure I own this data
        dir_map[dir].push_back(boxkey);
    }
}

void LevelPartitions::FormMaps() {
#ifndef RELEASE
    CallStackEntry entry("LevelPartitions::FormMaps");
#endif
    // Outgoing (upwards pass)
    for (int i = 0; i < static_cast<int>(_hf_vecs_out.size()); ++i) {
        InsertIntoDirMap(_level_hdmap_out[i], _hf_vecs_out[i].lclmap());
    }
    InsertIntoDirMap(_level_hdmap_out[unit_level_], _unit_vec.lclmap());
    // Incoming (downwards pass)
    for (int i = 0; i < static_cast<int>(_hf_vecs_inc.size()); ++i) {
        InsertIntoDirMap(_level_hdmap_inc[i], _hf_vecs_inc[i].lclmap());
    }
    InsertIntoDirMap(_level_hdmap_inc[unit_level_], _unit_vec.lclmap());
    // Remove old data.
    for (int i = 0; i < static_cast<int>(_hdkeys_out.size()); ++i) {
        std::vector<BoxAndDirKey>& old_keys = _hdkeys_out[i];
        std::vector<BoxAndDirKey>().swap(old_keys);
    }
    for (int i = 0; i < static_cast<int>(_hdkeys_inc.size()); ++i) {
        std::vector<BoxAndDirKey>& old_keys = _hdkeys_inc[i];
        std::vector<BoxAndDirKey>().swap(old_keys);
    }
}

BoxAndDirDat& LevelPartitions::Access(BoxAndDirKey key, bool out) {
#ifndef RELEASE
    CallStackEntry entry("LevelPartitions::Access");
#endif
    int level = key._boxkey._level;
    if (level == unit_level_) {
        return _unit_vec.access(key);
    }
    if (out) {
        return _hf_vecs_out[level].access(key);
    }
    return _hf_vecs_inc[level].access(key);
}

std::pair<bool, BoxAndDirDat&> LevelPartitions::SafeAccess(BoxAndDirKey key,
                                                           bool out) {
#ifndef RELEASE
    CallStackEntry entry("LevelPartitions::SafeAccess");
#endif
    int level = key._boxkey._level;
    if (level == unit_level_) {
        return _unit_vec.contains(key);
    }
    if (out) {
        return _hf_vecs_out[level].contains(key);
    }
    return _hf_vecs_inc[level].contains(key);
}

int LevelPartitions::Owner(BoxAndDirKey key, bool out) {
#ifndef RELEASE
    CallStackEntry entry("LevelPartitions::Owner");
#endif
    int level = key._boxkey._level;
    if (level == unit_level_) {
        return _unit_vec.prtn().owner(key);
    }
    if (out) {
        return _hf_vecs_out[level].prtn().owner(key);
    }
    return _hf_vecs_inc[level].prtn().owner(key);
}

int LevelPartitions::Owner(BoxKey key) {
#ifndef RELEASE
    CallStackEntry entry("LevelPartitions::Owner");
#endif
    return _lf_boxvec.prtn().owner(key);
}

void CleanDirVecMap(std::map<Index3, std::vector<BoxKey> >& curr_map) {
#ifndef RELEASE
    CallStackEntry entry("CleanDirVecMap");
#endif
    for (std::map<Index3, std::vector<BoxKey> >::iterator mi = curr_map.begin();
         mi != curr_map.end(); ++mi) {
        std::vector<BoxKey>().swap(mi->second);
    }
    curr_map.clear();
}

void CleanBoxAndDirMap(std::map<BoxAndDirKey, BoxAndDirDat>& curr_map) {
#ifndef RELEASE
    CallStackEntry entry("CleanBoxAndDirMap");
#endif
    for (std::map<BoxAndDirKey, BoxAndDirDat>::iterator mi = curr_map.begin();
         mi != curr_map.end(); ++mi) {
        mi->second._dirupeqnden.resize(0);
        mi->second._dirdnchkval.resize(0);
        std::vector<BoxAndDirKey>().swap(mi->second.interactionlist());
    }
    curr_map.clear();
}

int Wave3d::CleanLevel(int level) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::CleanLevel");
#endif
    CleanDirVecMap(_level_prtns._level_hdmap_inc[level]);
    CleanDirVecMap(_level_prtns._level_hdmap_out[level]);
    if (level == UnitLevel()) {
        CleanBoxAndDirMap(_level_prtns._unit_vec.lclmap());
    } else {
        CleanBoxAndDirMap(_level_prtns._hf_vecs_out[level].lclmap());
        CleanBoxAndDirMap(_level_prtns._hf_vecs_inc[level].lclmap());
    }
    return 0;
}

int Wave3d::CleanBoxvec() {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::CleanBoxvec");
#endif
    for (std::map<BoxKey, BoxDat>::iterator mi = _boxvec.lclmap().begin();
         mi != _boxvec.lclmap().end(); ++mi) {
        BoxDat& dat = mi->second;
        std::vector<BoxKey>().swap(dat._undeidxvec);
        std::vector<BoxKey>().swap(dat._vndeidxvec);
        std::vector<BoxKey>().swap(dat._wndeidxvec);
        std::vector<BoxKey>().swap(dat._xndeidxvec);
        std::vector<BoxKey>().swap(dat._endeidxvec);
        for (std::map<Index3, std::vector<BoxKey> >::iterator mi2 = dat._fndeidxvec.begin();
           mi2 != dat._fndeidxvec.end(); ++mi2) {
            std::vector<BoxKey>().swap(mi2->second);
        }
        dat._fndeidxvec.clear();
    }
    _boxvec.lclmap().clear();
    return 0;
}

template <typename T>
int FindInd(std::vector<T>& partition, std::vector<T>& end_partition, T key) {
    CHECK_TRUE_MSG(partition.size() != 0, "Missing partition!");
    int ind = std::lower_bound(partition.begin(),
                               partition.end(), key) - partition.begin();
    --ind;
    if (ind < static_cast<int>(partition.size()) - 1 &&
        key == partition[ind + 1]) {
        ++ind;
    }
    if (ind == partition.size() || ind == -1 ||
        key > end_partition[ind]) {
        return -1;
    }
    return ind;
}

int BoxAndDirLevelPrtn::owner(BoxAndDirKey& key) {
#ifndef RELEASE
    CallStackEntry entry("BoxAndDirLevelPrtn::owner");
#endif
    return FindInd(partition_, end_partition_, key);
}

int UnitLevelBoxPrtn::owner(BoxAndDirKey& key) {
#ifndef RELEASE
    CallStackEntry entry("UnitLevelBoxPrtn::owner");
#endif
    BoxKey boxkey = key._boxkey;
    return FindInd(partition_, end_partition_, boxkey);
}

int LowFreqBoxPrtn::owner(BoxKey& key) {
#ifndef RELEASE
    CallStackEntry entry("LowFreqBoxPrtn::owner");
#endif
    int level = key._level;
    CHECK_TRUE_MSG(level >= unit_level_, "bad level");
    int factor = pow2(level - unit_level_);
    Index3 idx = key._index;
    idx = idx / factor;
    BoxKey new_key(unit_level_, idx);
    return FindInd(partition_, end_partition_, new_key);
}
