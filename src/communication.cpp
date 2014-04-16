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
#include "wave3d.hpp"

#include <set>
#include <vector>

int Wave3d::HighFreqInteractionListKeys(Index3 dir, std::vector<BoxKey>& target_boxes,
                                        std::set<BoxAndDirKey>& reqbndset) {
#ifndef RELEASE
  CallStackEntry entry("Wave3d::HighFreqInteractionListKeys");
#endif
    for (int k = 0; k < target_boxes.size(); ++k) {
        BoxAndDirKey trgkey = BoxAndDirKey(target_boxes[k], dir);
        BoxAndDirDat& trgdat = _bndvec.access(trgkey);
        std::vector<BoxAndDirKey>& interaction_list = trgdat.interactionlist();
        for (int i = 0; i < interaction_list.size(); ++i) {
            reqbndset.insert(interaction_list[i]);
        }
    }
    return 0;
}

int Wave3d::HighFreqInteractionListKeys(int level,
                                        std::set<BoxAndDirKey>& request_keys) {
#ifndef RELEASE
  CallStackEntry entry("Wave3d::HighFreqInteractionListKeys");
#endif
    if (level == UnitLevel()) {
        // Loop over all boxes.  Only incoming computations have a non-empty
        // interaction list.
        std::map<BoxAndDirKey, BoxAndDirDat>& lclmap = _level_prtns._unit_vec.lclmap();
        for (std::map<BoxAndDirKey, BoxAndDirDat>::iterator mi = lclmap.begin();
             mi != lclmap.end(); ++mi) {
            std::vector<BoxAndDirKey>& interaction_list = mi->second.interactionlist();
            for (int j = 0; j < interaction_list.size(); ++j) {
                request_keys.insert(interaction_list[j]);
            }
        }
    } else {
        LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_inc[level];
        std::vector<BoxAndDirKey>& keys_inc = _level_prtns._hdkeys_inc[level];
        for (int i = 0; i < keys_inc.size(); ++i) {
            BoxAndDirKey key = keys_inc[i];
            BoxAndDirDat dat = vec.access(key);
            std::vector<BoxAndDirKey>& interaction_list = dat.interactionlist();
            for (int j = 0; j < interaction_list.size(); ++j) {
                request_keys.insert(interaction_list[j]);
            }
        }
    }
    return 0;
}

int Wave3d::LowFreqDownwardComm(std::set<BoxKey>& reqboxset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqDownwardComm");
#endif
    time_t t0 = time(0);
    std::vector<BoxKey> reqbox;
    reqbox.insert(reqbox.begin(), reqboxset.begin(), reqboxset.end());
    std::vector<int> mask(BoxDat_Number,0);
    mask[BoxDat_extden] = 1;
    mask[BoxDat_upeqnden] = 1;
    _boxvec.initialize_data();
    SAFE_FUNC_EVAL( _boxvec.getBegin(reqbox, mask) );
    SAFE_FUNC_EVAL( _boxvec.getEnd(mask) );
    time_t t1 = time(0);
    PrintParData(GatherParData(t0, t1), "Low frequency downward communication");
    PrintCommData(GatherCommData(_boxvec.kbytes_received()),
                  "kbytes received");
    PrintCommData(GatherCommData(_boxvec.kbytes_sent()),
                  "kbytes sent");
    return 0;
}

int Wave3d::GatherDensities(std::vector<int>& reqpts,
                            ParVec<int, cpx, PtPrtn>& den) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::GatherDensities");
#endif
    int mpirank = getMPIRank();
    ParVec<int, Point3, PtPrtn>& pos = (*_posptr);
    // Go through posptr to get nonlocal points
    for(std::map<int, Point3>::iterator mi = pos.lclmap().begin();
        mi != pos.lclmap().end(); ++mi) {
        reqpts.push_back( mi->first );
    }
    std::vector<int> all(1, 1);
    time_t t0 = time(0);
    SAFE_FUNC_EVAL( den.getBegin(reqpts, all) );
    SAFE_FUNC_EVAL( den.getEnd(all) );
    time_t t1 = time(0);
    if (mpirank == 0) {
        std::cout << "Density communication: " << difftime(t1, t0)
                  << " secs" << std::endl;
    }
    return 0;
}

int Wave3d::GatherDensities2(ParVec<int, cpx, PtPrtn>& den) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::GatherDensities");
#endif
    int mpirank = getMPIRank();

    std::set<int> req_dens;
    for (std::map<BoxKey,BoxDat>::iterator mi = _level_prtns._lf_boxvec.lclmap().begin();
        mi != _level_prtns._lf_boxvec.lclmap().end(); ++mi) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (HasPoints(curdat) &&
            _level_prtns._lf_boxvec.prtn().owner(curkey) == mpirank &&
            IsLeaf(curdat)) {
            std::vector<int>& curpis = curdat.ptidxvec();
            for (int k = 0; k < curpis.size(); ++k) {
                int poff = curpis[k];
                req_dens.insert(poff);
            }
        }
    }

    // Go through posptr to get nonlocal points
    std::vector<int> reqpts;
    reqpts.insert(reqpts.begin(), req_dens.begin(), req_dens.end());
    std::vector<int> all(1, 1);
    time_t t0 = time(0);
    SAFE_FUNC_EVAL( den.getBegin(reqpts, all) );
    SAFE_FUNC_EVAL( den.getEnd(all) );
    time_t t1 = time(0);
    if (mpirank == 0) {
        std::cout << "Density communication: " << difftime(t1, t0)
                  << " secs" << std::endl;
    }
    return 0;
}

int Wave3d::AllChildrenKeys(LevelBoxAndDirVec& vec,
                            std::vector<BoxAndDirKey>& req_keys) {
#ifndef RELEASE
    CallStackEntry entry("AllChildrenKeys");
#endif
    int mpirank = getMPIRank();
    for (std::map<BoxAndDirKey, BoxAndDirDat>::iterator mi = vec.lclmap().begin();
        mi != vec.lclmap().end(); ++mi) {
        BoxAndDirKey key = mi->first;
        if (vec.prtn().owner(key) == mpirank) {
            Index3 dir = key._dir;
            BoxKey boxkey = key._boxkey;
            Index3 pdir = ParentDir(dir);
            for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                int a = CHILD_IND1(ind);
                int b = CHILD_IND2(ind);
                int c = CHILD_IND3(ind);             
                BoxKey child_boxkey = ChildKey(boxkey, Index3(a, b, c));
                // We need the child key.
                req_keys.push_back(BoxAndDirKey(child_boxkey, pdir));
            }
        }
    }
    return 0;
}

int Wave3d::HighFreqL2LLevelCommPre(int level) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqL2LLevelCommPre");
#endif
    int childlevel = level + 1;
    // If the child level is at or above the starting level, then there is
    // no communication to do.  Also, if the child level is in the
    // low-frequency regime, there is no work to do.
    if (childlevel <= _starting_level || childlevel > UnitLevel()) {
        return 0;
    }

    std::vector<int> mask(BoxAndDirDat_Number, 0);
    mask[BoxAndDirDat_dirupeqnden] = 1;
    if (childlevel == UnitLevel()) {
        ParVec<BoxAndDirKey, BoxAndDirDat, UnitLevelBoxPrtn>& vec = _level_prtns._unit_vec;
        SAFE_FUNC_EVAL(vec.getBegin(&Wave3d::HighFreqL2LDataUp_wrapper, mask));
        SAFE_FUNC_EVAL(vec.getEnd(mask));
    } else {
        LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_inc[childlevel];
        SAFE_FUNC_EVAL(vec.getBegin(&Wave3d::HighFreqL2LDataUp_wrapper, mask));
        SAFE_FUNC_EVAL(vec.getEnd(mask));
    }
    return 0;
}

int Wave3d::HighFreqL2LDataUp(BoxAndDirKey key, BoxAndDirDat& dat,
                              std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqL2LDataUp");
#endif
    BoxKey parkey = ParentKey(key._boxkey);
    int parlevel = parkey.first;
    // Child directions are directions associated with the parent box.
    std::vector<Index3> dirs = ChildDir(key._dir);
    for (int i = 0; i < dirs.size(); ++i) {
        BoxAndDirKey new_key(parkey, dirs[i]);
        int pid = _level_prtns._hf_vecs_inc[parlevel].prtn().owner(new_key);
        if (0 <= pid && pid < getMPISize()) {
            pids.push_back(pid);
        }
    }
    return 0;
}

int Wave3d::HighFreqL2LDataUp_wrapper(BoxAndDirKey key, BoxAndDirDat& dat,
                                      std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqL2LDataUp_wrapper");
#endif
    return (Wave3d::_self)->HighFreqL2LDataUp(key, dat, pids);
}

int Wave3d::HighFreqL2LLevelCommPost(int level,
				     std::vector<BoxAndDirKey>& keys_affected) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqL2LLevelCommPost");
#endif
    // If the child level is at or above the starting level, then there is
    // no communication to do.  Also, if the child level is in the
    // low-frequency regime, there is no work to do.
    int childlevel = level + 1;
    if (childlevel <= _starting_level || childlevel > UnitLevel()) {
        return 0;
    }
    CHECK_TRUE(level <= UnitLevel());
    std::vector<int> mask(BoxAndDirDat_Number, 0);
    mask[BoxAndDirDat_dirdnchkval] = 1;
    if (childlevel == UnitLevel()) {
        ParVec<BoxAndDirKey, BoxAndDirDat, UnitLevelBoxPrtn>& vec = _level_prtns._unit_vec;
	// Send data that I have
	vec.putBegin(keys_affected, mask);
	vec.putEnd(mask);
    } else {
        LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_inc[childlevel];
	// Send data that I have
	vec.putBegin(keys_affected, mask);
	vec.putEnd(mask);
    }
    return 0;
}

int Wave3d::HighFreqM2LComm(std::set<BoxAndDirKey>& reqbndset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2LComm");
#endif
    std::vector<int> mask(BoxAndDirDat_Number, 0);
    mask[BoxAndDirDat_dirupeqnden] = 1;
    std::vector<BoxAndDirKey> reqbnd;
    reqbnd.insert(reqbnd.begin(), reqbndset.begin(), reqbndset.end());
    SAFE_FUNC_EVAL( _bndvec.getBegin(reqbnd, mask) );
    SAFE_FUNC_EVAL( _bndvec.getEnd(mask) );
    return 0;
}

int Wave3d::HighFreqM2LComm(int level,
                            std::set<BoxAndDirKey>& request_keys) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2LComm");
#endif
    std::vector<int> mask(BoxAndDirDat_Number, 0);
    mask[BoxAndDirDat_dirupeqnden] = 1;
    std::vector<BoxAndDirKey> req;
    req.insert(req.begin(), request_keys.begin(), request_keys.end());
    if (level == UnitLevel()) {
        ParVec<BoxAndDirKey, BoxAndDirDat, UnitLevelBoxPrtn>& vec = _level_prtns._unit_vec;
        SAFE_FUNC_EVAL( vec.getBegin(req, mask) );
        SAFE_FUNC_EVAL( vec.getEnd(mask) );
    } else {
        LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_out[level];
        SAFE_FUNC_EVAL( vec.getBegin(req, mask) );
        SAFE_FUNC_EVAL( vec.getEnd(mask) );
    }
    return 0;
}

int Wave3d::HighFreqM2MLevelComm(int level) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2MLevelComm");
#endif
    // If we are at or above the starting level, there is no communication to
    // do.  (There is no computation at the level above).
    if (level >= _starting_level) {
        return 0;
    }
    std::vector<int> mask(BoxAndDirDat_Number, 0);
    mask[BoxAndDirDat_dirupeqnden] = 1;
    if (level == UnitLevel()) {
        ParVec<BoxAndDirKey, BoxAndDirDat, UnitLevelBoxPrtn>& vec = _level_prtns._unit_vec;
        SAFE_FUNC_EVAL(vec.getBegin(&Wave3d::HighFreqM2MDataUp_wrapper, mask));
        SAFE_FUNC_EVAL(vec.getEnd(mask));
    } else {
        LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_out[level];
        SAFE_FUNC_EVAL(vec.getBegin(&Wave3d::HighFreqM2MDataUp_wrapper, mask));
        SAFE_FUNC_EVAL(vec.getEnd(mask));
    }
    return 0;
}

int Wave3d::HighFreqM2MDataUp(BoxAndDirKey key, BoxAndDirDat& dat,
                              std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2MDataUp");
#endif
    BoxKey parkey = ParentKey(key._boxkey);
    int parlevel = parkey.first;
    // Child directions are directions associated with the parent box.
    std::vector<Index3> dirs = ChildDir(key._dir);
    for (int i = 0; i < dirs.size(); ++i) {
        BoxAndDirKey new_key(parkey, dirs[i]);
        int pid = _level_prtns._hf_vecs_out[parlevel].prtn().owner(new_key);
        if (0 <= pid && pid < getMPISize()) {
            pids.push_back(pid);
        }
    }
    return 0;
}

int Wave3d::HighFreqM2MDataUp_wrapper(BoxAndDirKey key, BoxAndDirDat& dat,
                                      std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2MDataUp_wrapper");
#endif
    return (Wave3d::_self)->HighFreqM2MDataUp(key, dat, pids);
}
