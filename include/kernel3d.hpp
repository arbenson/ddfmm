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
#ifndef _KERNEL3D_HPP_
#define _KERNEL3D_HPP_

#include "nummat.hpp"
#include "vec3t.hpp"

enum {
    KERNEL_HELM_SINGLE = 0,
    KERNEL_HELM_DOUBLE = 1,
    KERNEL_HELM_MIXED = 2,
    KERNEL_EXPR = 3
};

class Kernel3d {
protected:
  static double _mindif;
  int _type;
public:
    Kernel3d(int t=KERNEL_HELM_SINGLE): _type(t) {;}
    ~Kernel3d() {;}
    int& type() { return _type; }
    
    int kernel(const DblNumMat& trgpos, const DblNumMat& srcpos,
	       const DblNumMat& srcnor, CpxNumMat& mat);
    
    bool NeedsNormals () { return _type == KERNEL_HELM_MIXED ||
	                          _type == KERNEL_HELM_DOUBLE; }
};

#endif
