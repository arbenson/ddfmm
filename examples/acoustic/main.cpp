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

int main(int argc, char** argv) {
  std::cout << "starting!" << std::endl;
  //
  srand48(time(NULL));
  map<string,string> opts;
  optionsCreate(argc, argv, opts);
  //1. read stuff
  map<string,string>::iterator mi;
  vector<int> all(1,1);
  //
  mi = opts.find("-vertfile");  assert(mi!=opts.end());
  char vertfile[100];
  {
    istringstream ss((*mi).second);
    ss>>vertfile;
  }
  vector<Point3> vertvec;
  {
    ifstream fin(vertfile);
    SAFE_FUNC_EVAL( deserialize(vertvec, fin, all) );
  }
  //
  mi = opts.find("-facefile");
  CHECK_TRUE(mi != opts.end());
  char facefile[100];
  {
    istringstream ss((*mi).second);
    ss>>facefile;
  }
  vector<Index3> facevec;
  {
    ifstream fin(facefile);
    SAFE_FUNC_EVAL( deserialize(facevec, fin, all) );
  }
  //
  int N = vertvec.size();
  //
  Point3 ctr(0,0,0);
  //
  mi = opts.find("-accu");
  assert(mi != opts.end());
  int accu;
  {
    istringstream ss((*mi).second);
    ss>>accu;
  }
  CHECK_TRUE( accu>=1 && accu<=3 );
  mi = opts.find("-kt");
  CHECK_TRUE(mi != opts.end());
  int kt;
  {
    istringstream ss((*mi).second);
    ss>>kt;
  } 
  //CHECK_TRUE( kt==KNL_HELM_DOUB || kt==KNL_HELM_COMB);
  Kernel3d knlbie(kt);

  std::cout << "facevec size: " << facevec.size() << std::endl;
  std::cout << "vertvec size: " << vertvec.size() << std::endl;
  
  //2. scattering
  Acoustic3D acou("acou3d");
  std::cout << "setting up..." << std::endl;
  SAFE_FUNC_EVAL( acou.setup(vertvec, facevec, ctr, accu, knlbie) );
  std::cout << "done setting up..." << std::endl;
  //
  //load den
  mi = opts.find("-denfile");  assert(mi!=opts.end());
  char denfile[100];
  {
    istringstream ss((*mi).second);
    ss>>denfile;
  }
  vector<cpx> denvec;
  {
    ifstream fin(denfile);
    SAFE_FUNC_EVAL( deserialize(denvec, fin, all) );
  }
  //
  mi = opts.find("-chkfile");  assert(mi!=opts.end());
  char chkfile[100];
  {
    istringstream ss((*mi).second);
    ss>>chkfile;
  }
  vector<Point3> chkvec;
  {
    ifstream fin(chkfile);
    SAFE_FUNC_EVAL( deserialize(chkvec, fin, all) );
  }
  int C = chkvec.size();
  //
  vector<cpx> valvec(C, cpx(0,0));
  std::cout << "evaling..." << std::endl;
  SAFE_FUNC_EVAL( acou.eval(chkvec, denvec, valvec) );
  std::cout << "done evaling..." << std::endl;
  //
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
  //
  return 0;
}
