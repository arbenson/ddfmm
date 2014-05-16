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
  MPI_Init(&argc, &argv);
  int mpirank, mpisize;
  getMPIInfo(&mpirank, &mpisize);

  if (mpirank == 0) {
    std::cout << "starting!" << std::endl;
  }
  //
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
  {
    ifstream fin(vertfile);
    if (fin.fail()) {
      std::cerr << "Failed to open vertfile!" << std::endl;
      return -1;
    }
    SAFE_FUNC_EVAL( deserialize(vertvec, fin, all) );
  }

  std::string facefile = findOption(opts, "-facefile");
  if (facefile.empty()) {
    return -1;
  }
  vector<Index3> facevec;
  {
    ifstream fin(facefile);
    if (fin.fail()) {
      std::cerr << "Failed to open vertfile!" << std::endl;
      return -1;
    }
    SAFE_FUNC_EVAL( deserialize(facevec, fin, all) );
  }

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
  {
    istringstream ss(opt);
    ss >> accuracy;
  }
  CHECK_TRUE(accuracy >= 1 && accuracy <= 3);

  //2. scattering
  Acoustic3d acou;
  SAFE_FUNC_EVAL( acou.setup(vertvec, facevec, ctr, accuracy) );
  //
  //load den
  mi = opts.find("-denfile");  assert(mi!=opts.end());
  char denfile[100];
  {
    istringstream ss((*mi).second);
    ss >> denfile;
  }
  vector<cpx> denvec;
  {
    ifstream fin(denfile);
    SAFE_FUNC_EVAL( deserialize(denvec, fin, all) );
  }

  mi = opts.find("-chkfile");  assert(mi!=opts.end());
  char chkfile[100];
  {
    istringstream ss((*mi).second);
    ss>>chkfile;
  }

  opt = findOption(opts, "-wave3d_K");
  if (opt.empty()) {
    return -1;
  }
  std::istringstream ss(opt);
  ss >> acou._K;

  vector<Point3> chkvec;
  vector<cpx> valvec;
#if 0
  {
    ifstream fin(chkfile);
    SAFE_FUNC_EVAL( deserialize(chkvec, fin, all) );
  }
  int C = chkvec.size();

  vector<cpx> valvec(C, cpx(0,0));
#endif
  if (mpirank == 0) {
    std::cout << "evaling..." << std::endl;
  }
  SAFE_FUNC_EVAL( acou.eval(chkvec, denvec, valvec, opts) );
  if (mpirank == 0) {
    std::cout << "done evaling..." << std::endl;
  }
  mi = opts.find("-valfile");
  assert(mi != opts.end());
  char valfile[100];
  {
    istringstream ss((*mi).second);
    ss>>valfile;
  }
  {
    ofstream fot(valfile);
    SAFE_FUNC_EVAL( serialize(valvec, fot, all) );
  }

  return 0;
}
