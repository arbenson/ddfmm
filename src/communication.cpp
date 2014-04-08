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

int Wave3d::HighFreqChildrenKeys(double W,
                                 std::map<Index3, std::vector<BoxKey> >& level_out,
				 std::vector<BoxAndDirKey>& children_keys) {
#ifndef RELEASE
  CallStackEntry entry("Wave3d::HighFreqChildrenKeys");
#endif
    double eps = 1e-12;
    if (abs(W - 1) < eps) {
        // We assume that all directions for a box at the unit level
        // are assigned to the same processor.  Thus, the processor owns
        // the subtree at that box, and no communication is needed for M2M.
        return 0;
    }

    for (std::map<Index3, std::vector<BoxKey> >::iterator mi = level_out.begin();
	 mi != level_out.end(); ++mi) {
        Index3 dir = mi->first;
	std::vector<BoxKey>& keys = mi->second;
        for (int i = 0; i < keys.size(); ++i) {
            // Pick the direction such that the child wedges in that direction
            // contain the parent wedge in direction dir
            BoxKey curr_key = keys[i];
	    Index3 pdir = ParentDir(dir);
	    for (int ind = 0; ind < NUM_CHILDREN; ++ind) {
                int a = CHILD_IND1(ind);
		int b = CHILD_IND2(ind);
		int c = CHILD_IND3(ind);
		BoxKey chdkey = ChildKey(curr_key, Index3(a, b, c));
		children_keys.push_back(BoxAndDirKey(chdkey, pdir));
	    }
	}
    }
    return 0;
}


int Wave3d::HighFreqInteractionListKeys(Index3 dir, std::vector<BoxKey>& target_boxes,
					std::set<BoxAndDirKey>& reqbndset) {
#ifndef RELEASE
  CallStackEntry entry("Wave3d::HighFreqInteractionListKeys");
#endif
  for (int k = 0; k < target_boxes.size(); ++k) {
      BoxKey trgkey = target_boxes[k];
      BoxDat& trgdat = _boxvec.access(trgkey);
      CHECK_TRUE(HasPoints(trgdat));
      std::vector<BoxKey>& tmpvec = trgdat.fndeidxvec()[dir];
      for (int i = 0; i < tmpvec.size(); ++i) {
          BoxKey srckey = tmpvec[i];
          reqbndset.insert(BoxAndDirKey(srckey, dir));
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
