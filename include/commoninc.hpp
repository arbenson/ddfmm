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
#ifndef _COMMONINC_HPP_
#define _COMMONINC_HPP_

//STL stuff
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>

#include <cfloat>
#include <cassert>
#include <cmath>
#include <string>
#include <complex>

#include <vector>
#include <set>
#include <map>
#include <deque>
#include <queue>
#include <utility>
#include <algorithm>

//external libraries
#include "fftw3.h"
#include "blas.h"

#include "mpi.h"

typedef std::complex<double> cpx;

//aux functions
inline int pow2(int l) { assert(l >= 0); return (1 << l); }

#define SAFE_FUNC_EVAL(fun)  { int ierr = fun; if (ierr != 0) { throw new std::exception(); } }
#define CHECK_TRUE(expr) { if((expr) == 0) { std::cerr << "wrong" << std::endl; throw new std::exception(); } }
#define CHECK_TRUE_MSG(expr, msg) { if((expr) == 0) { std::cerr << msg << std::endl; throw new std::exception(); } }
//#define CHECK_TRUE(expr) { if((expr) == 0) { assert(0); } }

template <class T, class S>
std::istream& operator>>(std::istream& is, std::pair<T, S>& a) {
  is >> a.first;
  is >> a.second;
  return is;
}
template <class T, class S>
std::ostream& operator<<(std::ostream& os, const std::pair<T, S>& a) {
  os << a.first << " " << a.second;
  return os;
}

int getMPIRank();
int getMPISize();
int getMPIInfo(int *mpirank, int *mpisize);

#ifndef RELEASE
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack( std::ostream& os=std::cerr );

class CallStackEntry {
public:
    CallStackEntry( std::string s ) { 
        if ( !std::uncaught_exception() )
            PushCallStack(s); 
    }
    ~CallStackEntry() { 
        if ( !std::uncaught_exception() ) 
            PopCallStack(); 
    }
};
#endif // ifndef RELEASE

#endif
