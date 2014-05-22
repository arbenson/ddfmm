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
#include "kernel3d.hpp"
#include "wave3d.hpp"

#include <set>
#include <vector>

int Wave3d::HighFreqInteractionListKeys(int level,
                                        std::set<BoxAndDirKey>& request_keys) {
#ifndef RELEASE
  CallStackEntry entry("Wave3d::HighFreqInteractionListKeys");
#endif
    std::map<BoxAndDirKey, BoxAndDirDat>& lclmap =
      (level == UnitLevel()) ? _level_prtns._unit_vec.lclmap()
                             : _level_prtns._hf_vecs_inc[level].lclmap();
    for (auto& kv : lclmap) {
        for (BoxAndDirKey& key : kv.second.interactionlist()) {
            request_keys.insert(key);
        }
    }
    return 0;
}

int Wave3d::LowFreqDownwardComm(std::set<BoxKey>& reqboxset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqDownwardComm");
#endif
    double t0 = MPI_Wtime();
    std::vector<BoxKey> reqbox;
    reqbox.insert(reqbox.begin(), reqboxset.begin(), reqboxset.end());
    std::vector<int> mask(BoxDat_Number, 0);
    mask[BoxDat_extden] = 1;
    mask[BoxDat_upeqnden] = 1;
    _level_prtns._lf_boxvec.initialize_data();
    SAFE_FUNC_EVAL( _level_prtns._lf_boxvec.getBegin(reqbox, mask) );
    SAFE_FUNC_EVAL( _level_prtns._lf_boxvec.getEnd(mask) );
    double t1 = MPI_Wtime();
    PrintParData(GatherParData(t0, t1), "Low frequency downward communication");
    PrintCommData(GatherCommData( _level_prtns._lf_boxvec.kbytes_sent()),
                  "kbytes sent");
    return 0;
}

int Wave3d::GatherDensities(ParVec<int, cpx, PtPrtn>& den) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::GatherDensities");
#endif
    int mpirank = getMPIRank();

    std::set<int> req_dens;
    for (auto& kv : _level_prtns._lf_boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        if (HasPoints(curdat) &&
            _level_prtns.Owner(curkey) == mpirank &&
            IsLeaf(curdat)) {
            std::vector<int>& curpis = curdat.ptidxvec();
            for (int k = 0; k < static_cast<int>(curpis.size()); ++k) {
                int poff = curpis[k];
                req_dens.insert(poff);
            }
        }
    }

    // Go through posptr to get nonlocal points
    std::vector<int> reqpts;
    reqpts.insert(reqpts.begin(), req_dens.begin(), req_dens.end());
    std::vector<int> all(1, 1);
    double t0 = MPI_Wtime();
    SAFE_FUNC_EVAL( den.getBegin(reqpts, all) );
    SAFE_FUNC_EVAL( den.getEnd(all) );
    double t1 = MPI_Wtime();
    if (mpirank == 0) {
#ifdef _VERBOSE_
        std::cout << "Density communication: " << MPIDiffTime(t0, t1)
                  << " secs" << std::endl;
#endif
    }
    if (_kernel.type() == KERNEL_HELM_MIXED) {
        // We also need to communicate the normal vectors.
        double t0 = MPI_Wtime();
	SAFE_FUNC_EVAL(_normal_vecs.getBegin(reqpts, all));
	SAFE_FUNC_EVAL(_normal_vecs.getEnd(all));
	double t1 = MPI_Wtime();
	if (mpirank == 0) {
#ifdef _VERBOSE_
	    std::cout << "Normal vector communication: " << MPIDiffTime(t0, t1)
		      << " secs" << std::endl;
#endif
	}
	
    }

    // Now set the densities.
    for (auto& kv : _level_prtns._lf_boxvec.lclmap()) {
        BoxKey curkey = kv.first;
        BoxDat& curdat = kv.second;
        if (HasPoints(curdat) &&
            _level_prtns.Owner(curkey) == mpirank &&
            IsLeaf(curdat)) {
            std::vector<int>& curpis = curdat.ptidxvec();
            CpxNumVec& extden = curdat.extden();
            extden.resize(curpis.size());
	    if (_kernel.type() == KERNEL_HELM_MIXED) {
                // Make room for the normal vectors
	        DblNumMat& extnor = curdat.extnor();
		extnor.resize(3, curpis.size());
	    }
            for (int k = 0; k < static_cast<int>(curpis.size()); ++k) {
                int poff = curpis[k];
                extden(k) = den.access(poff);
		if (_kernel.type() == KERNEL_HELM_MIXED) {
                    // Put in normal vectors.
                    Point3 nor = _normal_vecs.access(poff);
		    for (int d = 0; d < 3; ++d) {
                        curdat.extnor()(d, k) = nor(d);
		    }
		}
            }
        }
    }

    // Discard the points now that they have been converted to
    // equivalent densities.
    den.discard(reqpts);
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

    // The mask is zero because we are just getting the necessary keys
    std::vector<int> mask(BoxAndDirDat_Number, 0);
    auto func = [this] (BoxAndDirKey key, BoxAndDirDat& dat,
			std::vector<int>& pids) {
	return HighFreqL2LDataUp(key, dat, pids);
    };
    if (childlevel == UnitLevel()) {
        ParVec<BoxAndDirKey, BoxAndDirDat, UnitLevelBoxPrtn>& vec = _level_prtns._unit_vec;
	vec.initialize_data();
        SAFE_FUNC_EVAL(vec.getBegin(func, mask));
        SAFE_FUNC_EVAL(vec.getEnd(mask));
	PrintCommData(GatherCommData(vec.kbytes_sent()),
		      "L2L Pre: kbytes sent");
    } else {
        LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_inc[childlevel];
	vec.initialize_data();
        SAFE_FUNC_EVAL(vec.getBegin(func, mask));
        SAFE_FUNC_EVAL(vec.getEnd(mask));
	PrintCommData(GatherCommData(vec.kbytes_sent()),
		      "L2L Pre: kbytes sent");
    }
    return 0;
}

int Wave3d::HighFreqL2LLevelCommPost(int level,
                                     std::set<BoxAndDirKey>& affected_keys) {
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
    std::vector<BoxAndDirKey> send_keys;
    send_keys.insert(send_keys.begin(), affected_keys.begin(), affected_keys.end());
    if (childlevel == UnitLevel()) {
        ParVec<BoxAndDirKey, BoxAndDirDat, UnitLevelBoxPrtn>& vec = _level_prtns._unit_vec;
        // Send data that I have
	vec.initialize_data();
        vec.putBegin(send_keys, mask);
        vec.putEnd(mask);
	PrintCommData(GatherCommData(vec.kbytes_sent()),
		      "L2L Post: kbytes sent");
    } else {
        LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_inc[childlevel];
        // Send data that I have
	vec.initialize_data();
        vec.putBegin(send_keys, mask);
        vec.putEnd(mask);
	PrintCommData(GatherCommData(vec.kbytes_sent()),
		      "L2L Post: kbytes sent");
    }
    return 0;
}

int Wave3d::HighFreqL2LDataUp(BoxAndDirKey key, BoxAndDirDat& dat,
                              std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqL2LDataUp");
#endif
    BoxKey parkey = ParentKey(key._boxkey);
    // Child directions are directions associated with the parent box.
    std::vector<Index3> dirs = ChildDir(key._dir);
    for (int i = 0; i < static_cast<int>(dirs.size()); ++i) {
        BoxAndDirKey new_key(parkey, dirs[i]);
        int pid = _level_prtns.Owner(new_key, false);
        if (0 <= pid && pid < getMPISize()) {
            pids.push_back(pid);
        }
    }
    return 0;
}


int Wave3d::HighFreqM2LComm(int level,
                            std::set<BoxAndDirKey>& request_keys) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2LComm");
#endif
    if (level < _starting_level || level > UnitLevel()) {
        // No communication to do.
        return 0;
    }
    std::vector<int> mask(BoxAndDirDat_Number, 0);
    mask[BoxAndDirDat_dirupeqnden] = 1;
    std::vector<BoxAndDirKey> req;
    req.insert(req.begin(), request_keys.begin(), request_keys.end());
    if (level == UnitLevel()) {
        ParVec<BoxAndDirKey, BoxAndDirDat, UnitLevelBoxPrtn>& vec = _level_prtns._unit_vec;
	vec.initialize_data();
        SAFE_FUNC_EVAL( vec.getBegin(req, mask) );
        SAFE_FUNC_EVAL( vec.getEnd(mask) );
	PrintCommData(GatherCommData(vec.kbytes_sent()),
		      "M2L: kbytes sent");
    } else {
        LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_out[level];
	vec.initialize_data();
        SAFE_FUNC_EVAL( vec.getBegin(req, mask) );
        SAFE_FUNC_EVAL( vec.getEnd(mask) );
	PrintCommData(GatherCommData(vec.kbytes_sent()),
		      "M2L: kbytes sent");
    }
    return 0;
}

int Wave3d::HighFreqM2MLevelComm(int level) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2MLevelComm");
#endif
    // If we are at or above the starting level, there is no communication to
    // do.  (There is no computation at the level above).
    if (level <= _starting_level || level > UnitLevel()) {
        return 0;
    }
    std::vector<int> mask(BoxAndDirDat_Number, 0);
    mask[BoxAndDirDat_dirupeqnden] = 1;
    auto func = [this] (BoxAndDirKey key, BoxAndDirDat& dat,
			std::vector<int>& pids) {
	return HighFreqM2MDataUp(key, dat, pids);
    };
    if (level == UnitLevel()) {
        ParVec<BoxAndDirKey, BoxAndDirDat, UnitLevelBoxPrtn>& vec = _level_prtns._unit_vec;
	vec.initialize_data();	
        SAFE_FUNC_EVAL(vec.getBegin(func, mask));
        SAFE_FUNC_EVAL(vec.getEnd(mask));
	PrintCommData(GatherCommData(vec.kbytes_sent()),
		      "M2M: kbytes sent");
    } else {
        LevelBoxAndDirVec& vec = _level_prtns._hf_vecs_out[level];
	vec.initialize_data();
        SAFE_FUNC_EVAL(vec.getBegin(mask, func));
        SAFE_FUNC_EVAL(vec.getEnd(mask));
	PrintCommData(GatherCommData(vec.kbytes_sent()),
		      "M2M: kbytes sent");
    }
    return 0;
}

int Wave3d::HighFreqM2MDataUp(BoxAndDirKey key, BoxAndDirDat& dat,
                              std::vector<int>& pids) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqM2MDataUp");
#endif
    BoxKey parkey = ParentKey(key._boxkey);
    // Child directions are directions associated with the parent box.
    for (Index3 dir : ChildDir(key._dir)) {
        BoxAndDirKey new_key(parkey, dir);
        int pid = _level_prtns.Owner(new_key, true);
        if (0 <= pid && pid < getMPISize()) {
            pids.push_back(pid);
        }
    }
    return 0;
}
