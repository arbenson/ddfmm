/* Distributed Directional Fast Multipole Method
   Copyright (C) 2013 Austin Benson, Lexing Ying, and Jack Poulson

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
#ifndef _NUMVEC_HPP_
#define _NUMVEC_HPP_

#include "commoninc.hpp"

template <class F>
class NumVec
{
public:
    NumVec(int m=0): _m(m), _owndata(true)  {
#ifndef RELEASE
	CallStackEntry entry("NumVec::NumVec");
#endif
	allocate();
    }

    NumVec(int m, bool owndata, F* data): _m(m), _owndata(owndata) {
#ifndef RELEASE
	CallStackEntry entry("NumVec::NumVec");
#endif
	if (_owndata) {
	    allocate();
	    if (_m > 0) {
		for(int i=0; i<_m; i++) {
		    _data[i] = data[i];
		}
	    }
	} else {
	    _data = data;
	}
    }

    NumVec(const NumVec& C): _m(C._m), _owndata(C._owndata)  {
#ifndef RELEASE
	CallStackEntry entry("NumVec::NumVec");
#endif
	if (_owndata) {
	    allocate();
	    fill(C);
	} else {
	    _data = C._data;
	}
    }

    ~NumVec() {
#ifndef RELEASE
	CallStackEntry entry("NumVec::~NumVec");
#endif
	if (_owndata)
	    deallocate();
    }

    NumVec& operator=(const NumVec& C)  {
#ifndef RELEASE
	CallStackEntry entry("NumVec::operator=");
#endif
	if (_owndata)
	    deallocate();
	_m = C._m;
	_owndata=C._owndata;
	if (_owndata) {
	    allocate();
	    fill(C);
	} else {
	    _data =C._data;
	}
	return *this;
    }

    void resize(int m)  {
#ifndef RELEASE
	CallStackEntry entry("NumVec::resize");
#endif
	assert(_owndata);
	if(m !=_m) {
	    deallocate();
	    _m = m;
	    allocate();
	}
    }

    const F& operator()(int i) const  {
#ifndef RELEASE
	CallStackEntry entry("NumVec::operator()");
#endif
	assert(i>=0 && i<_m);
	return _data[i]; 
    }

    F& operator()(int i)  {
#ifndef RELEASE
	CallStackEntry entry("NumVec::operator()");
#endif
	assert(i>=0 && i<_m);
	return _data[i]; 
    }
  
    F* data() const { return _data; }
    int m () const { return _m; }

private:
    int  _m;
    bool _owndata;
    F* _data;

    void allocate() {
#ifndef RELEASE
	CallStackEntry entry("NumVec::allocate");
#endif
	if (_m > 0) {
	    _data = new F[_m];
	    assert(_data != NULL);
	} else {
	    _data = NULL;
	}
    }

    void deallocate() {
#ifndef RELEASE
	CallStackEntry entry("NumVec::deallocate");
#endif
	if (_m > 0) {
	    delete[] _data;
	    _data = NULL;
	}
    }

    void fill(const NumVec& C)  {
#ifndef RELEASE
	CallStackEntry entry("NumVec::fill");
#endif
	if (_m > 0) {
	    for(int i = 0; i < _m; i++) {
		_data[i] = C._data[i];
	    }
	}
    }
};

template <class F> inline std::ostream& operator<<( std::ostream& os,
                                                   const NumVec<F>& vec)
{
#ifndef RELEASE
    CallStackEntry entry("operator<<");
#endif
    os << vec.m() << std::endl;
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for(int i = 0; i < vec.m(); i++) {
        os << " " << vec(i);
    }
    os << std::endl;
    return os;
}

template <class F> inline void setvalue(NumVec<F>& vec, F val)
{
#ifndef RELEASE
    CallStackEntry entry("setvalue");
#endif
    for(int i=0; i<vec.m(); i++)
	vec(i) = val;
}
template <class F> inline double energy(NumVec<F>& vec)
{
#ifndef RELEASE
    CallStackEntry entry("energy");
#endif
    double sum = 0;
    for(int i=0; i<vec.m(); i++)
	sum += abs(vec(i)*vec(i));
    return sum;
}


typedef NumVec<bool>   BolNumVec;
typedef NumVec<int>    IntNumVec;
typedef NumVec<double> DblNumVec;
typedef NumVec<cpx>    CpxNumVec;


#endif


