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
#include "acou3d.hpp"
#include "serialize.hpp"

using namespace std;

int optionsCreate(int argc, char** argv, map<string,string>& options) {
  options.clear();
  for(int k = 1; k < argc; k = k + 2) {
    options[string(argv[k])] = string(argv[k+1]);
  }
  return 0;
}

std::string findOption(std::map<std::string, std::string>& opts,
		       std::string option) {
    std::map<std::string, std::string>::iterator mi = opts.find(option);
    if (mi == opts.end()) {
        std::cerr << "Missing option " << option << std::endl;
        return "";
    }
    return mi->second;
}

int main(int argc, char** argv) {
#ifndef RELEASE
    // When not in release mode, we catch all errors so that we can print the
    // manual call stack.
    try {
#endif
    MPI_Init(&argc, &argv);
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    srand48(time(NULL));
    map<string,string> opts;
    optionsCreate(argc, argv, opts);
    if (mpirank == 0) {
        for (auto& kv : opts) {
            std::cout << kv.first << ": " << kv.second << std::endl;
	}
    }

    map<string,string>::iterator mi;
    vector<int> all(1,1);
    std::string vertfile = findOption(opts, "-vertfile");
    if (vertfile.empty()) {
        return -1;
    }

    // Read the vertex file
    vector<Point3> vertvec;
    ifstream vert_input(vertfile);
    if (vert_input.fail()) {
        std::cerr << "Failed to open vertfile!" << std::endl;
        return -1;
    }
    SAFE_FUNC_EVAL( deserialize(vertvec, vert_input, all) );

    std::string facefile = findOption(opts, "-facefile");
    if (facefile.empty()) {
        return -1;
    }
    vector<Index3> facevec;
    ifstream face_input(facefile);
    if (face_input.fail()) {
        std::cerr << "Failed to open vertfile!" << std::endl;
        return -1;
    }
    SAFE_FUNC_EVAL( deserialize(facevec, face_input, all) );

    if (mpirank == 0) {
        std::cout << "facevec size: " << facevec.size() << std::endl;
	std::cout << "vertvec size: " << vertvec.size() << std::endl;
    }
  
    Point3 ctr(0, 0, 0);
    std::string opt = findOption(opts, "-wave3d_ACCU");
    if (opt.empty()) {
        return -1;
    }
    
    int accuracy = 0;
    istringstream acc_input(opt);
    acc_input >> accuracy;
    CHECK_TRUE(accuracy >= 1 && accuracy <= 3);

    Acoustic3d acou;
    SAFE_FUNC_EVAL( acou.setup(vertvec, facevec, ctr, accuracy) );

    opt = findOption(opts, "-wave3d_K");
    if (opt.empty()) {
        return -1;
    }
    std::istringstream K_input(opt);
    K_input >> acou._K;

    vector<Point3> chkvec;
    vector<cpx> valvec;
    if (mpirank == 0) {
        std::cout << "evaling..." << std::endl;
    }
    SAFE_FUNC_EVAL( acou.eval(valvec, opts) );
    if (mpirank == 0) {
        std::cout << "done evaling..." << std::endl;
    }
    std::string valfile = findOption(opts, "-valfile");
    if (valfile.empty()) {
        return -1;
    }
    ofstream output(valfile);
    SAFE_FUNC_EVAL( serialize(valvec, output, all) );
#ifndef RELEASE
    } catch( ... ) {
        int mpirank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
        std::cerr << "Process " << mpirank << " caught error." << std::endl;
        DumpCallStack();
    }
#endif
    MPI_Finalize();
    return 0;
}
