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
#ifndef _ACOU3D_HPP_
#define _ACOU3D_HPP_

#include "vec3t.hpp"
#include "numtns.hpp"
#include "kernel3d.hpp"
#include "wave3d.hpp"

#include <string>

class VertexData {
  Point3 pos;  // vertex position
  cpx den;     // density
};

//---------------------------------------------------------------------------
class Acoustic3d {
public:
    std::vector<Point3> _vertvec;
    std::vector<Index3> _facevec; // all triangle faces
    Point3 _ctr;
    int _accu;
    std::vector<double> _diavec; //diagonal (solid angle/4pi)
    std::vector<double> _arevec; //length of each segment
    Wave3d _wave;
    std::map<int, DblNumMat> _gauwgts;
    std::map<int, DblNumMat> _sigwgts;
    double _K;
    std::vector<int> _dist;
    std::vector<cpx> _boundary;
    std::vector<int> _vert_distrib;
    
    Acoustic3d() {}
    ~Acoustic3d() {}
    int setup(std::vector<Point3>& vertvec, std::vector<Index3>& facevec,
              Point3 ctr, int accu);
    int Apply(ParVec<int, cpx, PtPrtn>& in, ParVec<int, cpx, PtPrtn>& out);
    void Apply(CpxNumVec& x, CpxNumVec& y);
    void Run(std::map<std::string, std::string>& opts);

private:
    // Given an index, return true if and only if the process with rank
    // mpirank owns the index in the distributed vectors.
    bool Own(int index, int mpirank);

    // Initialize faces, vertices, normal vectors, wave structure.
    int InitializeData(std::map<std::string, std::string>& opts);

    // Handle diagonal component.  In the limit, this is 1/2, but since
    // we discretize, we need to account for the angle.
    void HandleAngle(ParVec<int, cpx, PtPrtn>& in, ParVec<int, cpx, PtPrtn>& out,
		     ParVec<int, cpx, PtPrtn>& potentials);

    // Remove nearby from the output.  We deal with these cases with the
    // singulairty correction.
    void RemoveNearby(ParVec<int, cpx, PtPrtn>& in, ParVec<int, cpx, PtPrtn>& out,
		      ParVec<int, cpx, PtPrtn>& densities);

    // Deal with integration over a corner singularity for each face.
    void SingularityCorrection(ParVec<int, cpx, PtPrtn>& in, ParVec<int, cpx, PtPrtn>& out);

    std::map<std::string, std::string> _opts;
};

#endif
