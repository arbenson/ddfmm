#ifndef _COMMONINC_HPP_
#define _COMMONINC_HPP_

//STL stuff
#include <iostream>
#include <fstream>
#include <sstream>

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

//complex number
using std::complex;
using std::pair;
using std::istream;
using std::ostream;

typedef complex<double> cpx;

//aux functions
inline int pow2(int l) { assert(l>=0); return (1<<l); }

#define iC(fun)  { int ierr=fun; assert(ierr==0); }
#define iA(expr) { if((expr)==0) { std::cerr<<"wrong"<<std::endl; assert(expr); } }

template <class T, class S>
istream& operator>>(istream& is, pair<T,S>& a)
{
  is>>a.first;
  is>>a.second;
  return is;
}
template <class T, class S>
ostream& operator<<(ostream& os, const pair<T,S>& a)
{
  os<<a.first<<" "<<a.second;
  return os;
}

int getMPIRank();
int getMPISize();
int getMPIInfo(int *mpirank, int *mpisize);

#ifndef RELEASE
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack( std::ostream& os=std::cerr );

class CallStackEntry 
{
public:
    CallStackEntry( std::string s ) 
    { 
        if( !std::uncaught_exception() )
            PushCallStack(s); 
    }
    ~CallStackEntry() 
    { 
        if( !std::uncaught_exception() )
            PopCallStack(); 
    }
};
#endif // ifndef RELEASE

#endif
