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
#include "trmesh.hpp"

#include <fstream>

int Acoustic3d::setup(vector<Point3>& vertvec, vector<Index3>& facevec,
                      Point3 ctr, int accu) {
#ifndef RELEASE
    CallStackEntry entry("Acoustic3d::setup");
#endif
    _vertvec = vertvec;
    _facevec = facevec;
    _ctr = ctr;
    _accu = accu;
    // Compute the diagonal scaling.
    TrMesh trmesh;
    SAFE_FUNC_EVAL(trmesh.setup(vertvec, facevec));
    SAFE_FUNC_EVAL(trmesh.compute_interior(_diavec));
    for (int k = 0; k < _diavec.size(); ++k) {
        _diavec[k] /= (4*M_PI);
    }
    SAFE_FUNC_EVAL( trmesh.compute_area(_arevec) );
    // Load the quadrature weights
    vector<int> all(1,1);
    std::ifstream gin("gauwgts.bin");
    CHECK_TRUE_MSG(!gin.fail(), "Could not open gauwgts.bin");
    SAFE_FUNC_EVAL( deserialize(_gauwgts, gin, all) );
    int mpirank = getMPIRank();
    if (mpirank == 0) {
        std::cerr << "gauwgts size " << _gauwgts.size() << std::endl;
    }
    std::ifstream lin("sigwgts.bin");
    CHECK_TRUE_MSG(!lin.fail(), "Could not open sigwgts.bin");
    SAFE_FUNC_EVAL( deserialize(_sigwgts, lin, all) );
    if (mpirank == 0) {
        std::cerr << "sigwgts size " << _sigwgts.size() << std::endl;
    }
    return 0;
}

bool Acoustic3d::Own(int index, int mpirank) {
    return index >= _dist[mpirank] && index < _dist[mpirank + 1];
}

int Acoustic3d::eval(vector<cpx>& val, std::map<std::string, std::string>& opts) {
#ifndef RELEASE
    CallStackEntry entry("Acoustic3d::eval");
#endif
  int mpirank, mpisize;
  getMPIInfo(&mpirank, &mpisize);

  CHECK_TRUE_MSG(_gauwgts.find(5) != _gauwgts.end(), "Problem with quadrature weights");
  DblNumMat& gauwgt = _gauwgts[5];
  int num_quad_points = gauwgt.m();
  std::vector<Point3> posvec;
  std::vector<Point3> norvec;
  std::vector<cpx> denvec;
  std::vector<cpx> valvec;

  // Determine distribution of face vectors.
  int total_points = _facevec.size() * num_quad_points + _vertvec.size();
  int avg_points = total_points / mpisize;
  int extra = total_points - avg_points * mpisize;
  _dist.resize(mpisize + 1);
  _dist[0] = 0;
  for (int i = 1; i < _dist.size(); ++i) {
      _dist[i] = _dist[i - 1] + avg_points;
      if (i - 1 < extra) {
          _dist[i] += 1;
      }
  }

  // USE ZERO DENSITY FOR NOW
  std::vector<cpx> den(_vertvec.size(), 0);

  for (int fi = 0; fi < _facevec.size(); ++fi) {
      // Only read if there is a chance this process will own.
      // TODO(arbenson): add a check here.
      Index3& face = _facevec[fi];
      // Get the three vertices of the face.
      Point3 pos0 = _vertvec[face(0)];
      Point3 pos1 = _vertvec[face(1)];
      Point3 pos2 = _vertvec[face(2)];

      double are = _arevec[fi];
      Point3 nor = cross(pos1 - pos0, pos2 - pos0);
      nor = nor / nor.l2();

      // Current density at the corners.
      cpx den0 = den[face(0)];
      cpx den1 = den[face(1)];
      cpx den2 = den[face(2)];

      for (int gi = 0; gi < num_quad_points; ++gi) {
          if (!Own(fi * num_quad_points + gi, mpirank)) {
	      continue;
	  }
          double loc0 = gauwgt(gi, 0);
          double loc1 = gauwgt(gi, 1);
          double loc2 = gauwgt(gi, 2);
          double wgt  = gauwgt(gi, 3);
          posvec.push_back(loc0 * pos0 + loc1 * pos1 + loc2 * pos2);
          norvec.push_back(nor);
          denvec.push_back((loc0 * den0 + loc1 * den1 + loc2 * den2) * (are * wgt));
      }
  }
  int num_faces = _facevec.size();
  for (int vi = 0; vi < _vertvec.size(); ++vi) {
      if (!Own(num_faces * num_quad_points + vi, mpirank)) {
          continue;
      }
      posvec.push_back(_vertvec[vi]);
      norvec.push_back(_vertvec[vi]);  // dummy
      denvec.push_back(cpx(0, 0));
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpirank == 0) {
      std::cout << "putting in ownerinfo" << std::endl;
  }

  // Positions, densities, potentials, and normals all follow this partitioning.
  ParVec<int, Point3, PtPrtn>& positions = _wave._positions;
  ParVec<int, Point3, PtPrtn>& normal_vecs = _wave._normal_vecs;
  positions.prtn().ownerinfo() = _dist;
  normal_vecs.prtn().ownerinfo() = _dist;

  ParVec<int, cpx, PtPrtn> densities;
  densities.prtn().ownerinfo() = _dist;
  ParVec<int, cpx, PtPrtn> potentials;
  potentials.prtn().ownerinfo() = _dist;

  int start_ind = _dist[mpirank];
  for (int i = 0; i < posvec.size(); ++i) {
      positions.insert(start_ind + i, posvec[i]);
  }
  for (int i = 0; i < norvec.size(); ++i) {
      normal_vecs.insert(start_ind + i, norvec[i]);
  }
  for (int i = 0; i < denvec.size(); ++i) {
      densities.insert(start_ind + i, denvec[i]);
  }
  cpx dummy_val(0, 0);
  for (int i = 0; i < denvec.size(); ++i) {
      potentials.insert(start_ind + i, dummy_val);
  }

  _wave._ctr = _ctr;
  _wave._ACCU = _accu;
  _wave._kernel = Kernel3d(KERNEL_HELM_MIXED);
  _wave._equiv_kernel = Kernel3d(KERNEL_HELM);

  // Deal with geometry partition.  For now, we just do a cyclic partition.
  // TODO(arbenson): be more clever about the partition.
  int num_levels = ceil(log(sqrt(_K)) / log(2));
  int num_cells = pow2(num_levels);
  _wave._geomprtn.resize(num_cells, num_cells, num_cells);
  int curr_proc = 0;
  for (int k = 0; k < num_cells; ++k) {
      for (int j = 0; j < num_cells; ++j) {
          for (int i = 0; i < num_cells; ++i) {
              _wave._geomprtn(i, j, k) = curr_proc;
              curr_proc = (curr_proc + 1) % mpisize;
          }
      }
  }

  Mlib3d& mlib = _wave._mlib;
  mlib._NPQ = _wave._NPQ;
  mlib._kernel = Kernel3d(KERNEL_HELM);
  mlib.setup(opts);

  std::cout << "Setting up the wave..." << std::endl;
  _wave.setup(opts);

  _wave.eval(densities, potentials);
  // TODO(arbenson): Put stuff into output vector.
  return 0;
}


#if 0
int Acou3d::ApplyMatrix(vector<cpx>& in, vector<cpx>& out) {
  CHECK_TRUE_MSG(_gauwgts.find(5) != _gauwgts.end(), "Problem with quadrature weights");
  DblNumMat& gauwgt = _gauwgts[5];
  CHECK_TRUE_MSG(_sigwgts.find(5) != _gauwgts.end(), "Problem with quadrature weights");
  DblNumMat& sigwgt = _sigwgts[5];

  int numgau = gauwgt.m();
  int numsig = sigwgt.m();
  int num_faces = _facevec.size();
  int num_vertices = _vertvec.size();

  // 1. Scale input with gaussian quadrature.
  std::vector<cpx> densities;
  for (int fi = 0; fi < _facevec.size(); fi++) {
      // Only handle if this face is owned by this process.
      if (!(fi >= _dist[mpirank] && fi < _dist[mpirank + 1])) {
          continue;
      }
      Index3& ind = _facevec[fi];
      double are = _arevec[fi];
      cpx den0 = in[ind(0)];
      cpx den1 = in[ind(1)];
      cpx den2 = in[ind(2)];
      for (int gi = 0; gi < numgau; ++gi) {
          double loc0 = gauwgt(gi,0);
	  double loc1 = gauwgt(gi,1);
	  double loc2 = gauwgt(gi,2);
	  double wgt  = gauwgt(gi,3);
	  densities.push_back((loc0*den0 + loc1*den1 + loc2*den2) * (are*wgt));
      }
  }
  iA( cnt==num_faces*numgau );

  // 2. Call directional FMM
  vector<cpx> aux(num_faces*numgau + num_vertices, cpx(0,0));
  time_t t0, t1;
  iC( _wave.eval(tmp,aux) );

  // 3. angle
  // TODO(arbenson): what is this doing?
  for (int vi = 0; vi < _vertvec.size(); ++vi) {
    ot[vi] = _diavec[vi]*in[vi] + aux[cnt];
    cnt++;
  }
  iA( cnt==num_faces*numgau+num_vertices );

  // 4. Remove nearby from the output
  // TODO(arbenson): what is this doing?
  for (int fi = 0; fi < _facevec.size(); fi++) {
      Index3& ind = _facevec[fi];
      DblNumMat srcpos(3, numgau, false, (double*)(&(_posvec[fi*numgau])));
      DblNumMat srcnor(3, numgau, false, (double*)(&(_norvec[fi*numgau])));
      vector<Point3> trgpostmp;
      trgpostmp.push_back( _vertvec[ ind(0) ] );
      trgpostmp.push_back( _vertvec[ ind(1) ] );
      trgpostmp.push_back( _vertvec[ ind(2) ] );
      DblNumMat trgpos(3,3,false,(double*)(&(trgpostmp[0])));
      CpxNumVec srcden(numgau,false,(cpx*)(&(tmp[fi*numgau])));
      CpxNumVec trgval(3);
      CpxNumMat mat;
      iC( _knlbie.kernel(trgpos, srcpos, srcnor, mat) );
      iC( zgemv(1.0, mat, srcden, 0.0, trgval) );
      
      ot[ ind(0) ] -= trgval(0);
      ot[ ind(1) ] -= trgval(1);
      ot[ ind(2) ] -= trgval(2);
  }

  // 5. Visit faces and add singularity correction.
  for (int fi = 0; fi < _facevec.size(); ++fi) {
      Index3& ind = _facevec[fi];
      Point3 pos0 = _vertvec[ind(0)];
      Point3 pos1 = _vertvec[ind(1)];
      Point3 pos2 = _vertvec[ind(2)];
      double are = _arevec[fi];
      Point3 nor = cross(pos1-pos0, pos2-pos0); 
      nor = nor / nor.l2();
      cpx den0 = in[ind(0)];
      cpx den1 = in[ind(1)];
      cpx den2 = in[ind(2)];
      //
      //0
      {
	vector<Point3> srcpostmp;
	vector<Point3> srcnortmp;
	vector<Point3> trgpostmp;
	vector<cpx> srcdentmp;
	for (int li = 0; li < numsig; ++li) {
	  double loc0 = sigwgt(li,0);
	  double loc1 = sigwgt(li,1);
	  double loc2 = sigwgt(li,2);
	  double wgt  = sigwgt(li,3);
	  Point3 pos = loc0*pos0 + loc1*pos1 + loc2*pos2;
	  cpx    den = (loc0*den0 + loc1*den1 + loc2*den2)*(are*wgt);
	  srcpostmp.push_back( pos );
	  srcnortmp.push_back( nor );
	  srcdentmp.push_back( den );
	}
	trgpostmp.push_back( pos0 );
	DblNumMat srcpos(3,numsig,false,(double*)(&(srcpostmp[0])));
	DblNumMat srcnor(3,numsig,false,(double*)(&(srcnortmp[0])));
	DblNumMat trgpos(3,1,false,(double*)(&(trgpostmp[0])));
	CpxNumVec srcden(numsig,false,(cpx*)(&(srcdentmp[0])));
	CpxNumVec trgval(1);
	CpxNumMat mat;
	iC( _knlbie.kernel(trgpos, srcpos, srcnor, mat) );
	iC( zgemv(1.0, mat, srcden, 0.0, trgval) );
	ot[ ind(0) ] += trgval(0);
      }
      //1
      {
	vector<Point3> srcpostmp;
	vector<Point3> srcnortmp;
	vector<Point3> trgpostmp;
	vector<cpx> srcdentmp;
	for(int li = 0; li < numsig; ++li) {
	  double loc0 = sigwgt(li,0);
	  double loc1 = sigwgt(li,1);
	  double loc2 = sigwgt(li,2);
	  double wgt  = sigwgt(li,3);
	  Point3 pos = loc0*pos1 + loc1*pos2 + loc2*pos0;
	  cpx    den = (loc0*den1 + loc1*den2 + loc2*den0)*(are*wgt);
	  srcpostmp.push_back( pos );
	  srcnortmp.push_back( nor );
	  srcdentmp.push_back( den );
	}
	trgpostmp.push_back( pos1 );
	DblNumMat srcpos(3,numsig,false,(double*)(&(srcpostmp[0])));
	DblNumMat srcnor(3,numsig,false,(double*)(&(srcnortmp[0])));
	DblNumMat trgpos(3,1,false,(double*)(&(trgpostmp[0])));
	CpxNumVec srcden(numsig,false,(cpx*)(&(srcdentmp[0])));
	CpxNumVec trgval(1);
	CpxNumMat mat;
	iC( _knlbie.kernel(trgpos, srcpos, srcnor, mat) );
	iC( zgemv(1.0, mat, srcden, 0.0, trgval) );
	ot[ ind(1) ] += trgval(0);
      }
      //2
      {
	vector<Point3> srcpostmp;
	vector<Point3> srcnortmp;
	vector<Point3> trgpostmp;
	vector<cpx> srcdentmp;
	for(int li = 0; li < numsig; ++li) {
	  double loc0 = sigwgt(li,0);
	  double loc1 = sigwgt(li,1);
	  double loc2 = sigwgt(li,2);
	  double wgt  = sigwgt(li,3);
	  Point3 pos = loc0*pos2 + loc1*pos0 + loc2*pos1;
	  cpx    den = (loc0*den2 + loc1*den0 + loc2*den1)*(are*wgt);
	  srcpostmp.push_back( pos );
	  srcnortmp.push_back( nor );
	  srcdentmp.push_back( den );
	}
	trgpostmp.push_back( pos2 );
	DblNumMat srcpos(3,numsig,false,(double*)(&(srcpostmp[0])));
	DblNumMat srcnor(3,numsig,false,(double*)(&(srcnortmp[0])));
	DblNumMat trgpos(3,1,false,(double*)(&(trgpostmp[0])));
	CpxNumVec srcden(numsig,false,(cpx*)(&(srcdentmp[0])));
	CpxNumVec trgval(1);
	CpxNumMat mat;
	iC( _knlbie.kernel(trgpos, srcpos, srcnor, mat) );
	iC( zgemv(1.0, mat, srcden, 0.0, trgval) );
	ot[ind(2)] += trgval(0);
      }
    }
    return 0;
}
#endif
