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
#ifndef _NUMMAT_HPP_
#define _NUMMAT_HPP_

#include "numvec.hpp"

template <class F>
class NumMat {
public:
    NumMat(int m=0, int n=0): _m(m), _n(n), _owndata(true) {
#ifndef RELEASE
        CallStackEntry entry("NumMat::NumMat");
#endif
	allocate();
    }

    NumMat(int m, int n, bool owndata, F* data): _m(m), _n(n),
                                                 _owndata(owndata) {
#ifndef RELEASE
        CallStackEntry entry("NumMat::NumMat");
#endif
        if (_owndata) {
	    allocate();
            if (ValidDimensions()) {
                for (int i = 0; i < _m * _n; ++i) {
                    _data[i] = data[i]; 
                }
            }
        } else {
            _data = data;
        }
    }

    NumMat(const NumMat& C): _m(C._m), _n(C._n), _owndata(C._owndata) {
#ifndef RELEASE
        CallStackEntry entry("NumMat::NumMat");
#endif
        if (_owndata) {
            allocate();
	    fill(C);
        } else {
            _data = C._data;
        }
    }

    ~NumMat() {
#ifndef RELEASE
        CallStackEntry entry("NumMat::~NumMat");
#endif
        if (_owndata)
	    deallocate();
    }

    NumMat& operator=(const NumMat& C) {
#ifndef RELEASE
        CallStackEntry entry("NumMat::operator=");
#endif
        if (_owndata)
	    deallocate();
        _m = C._m;
	_n = C._n;
	_owndata = C._owndata;
        if (_owndata) {
	    allocate();
	    fill(C);
        } else {
            _data = C._data;
        }
        return *this;
    }

    void resize(int m, int n)  {
#ifndef RELEASE
        CallStackEntry entry("NumMat::resize");
#endif
        assert(_owndata);
        if (_m != m || _n != n) {
            if (ValidDimensions()) {
                delete[] _data;
                _data = NULL;
            }
            _m = m;
            _n = n;
	    allocate();
        }
    }

    const F& operator()(int i, int j) const  { 
        assert( i >= 0 && i < _m && j >= 0 && j < _n );
        return _data[i + j * _m];
    }

    F& operator()(int i, int j)  { 
        assert( i >= 0 && i < _m && j >= 0 && j < _n );
        return _data[i + j * _m];
    }

    F* clmdata(int j) {
#ifndef RELEASE
        CallStackEntry entry("NumMat::clmdata");
#endif
        return &(_data[j*_m]);
    }

    F* data() const { return _data; }
    int m() const { return _m; }
    int n() const { return _n; }

private:
    int _m;
    int _n;
    bool _owndata;
    F* _data;

    bool ValidDimensions() const { return _m > 0 && _n > 0; }

    void allocate() {
#ifndef RELEASE
        CallStackEntry entry("NumMat::allocate");
#endif
	if (ValidDimensions()) {
	    _data = new F[_m * _n];
	    assert( _data != NULL );
	} else {
	    _data = NULL;
	}
    }

    void deallocate() {
#ifndef RELEASE
        CallStackEntry entry("NumMat::deallocate");
#endif
	if (ValidDimensions()) {
	    delete[] _data;
	    _data = NULL;
	}
    }

    void fill(const NumMat& C) {
#ifndef RELEASE
        CallStackEntry entry("NumMat::fill");
#endif
	if (ValidDimensions()) {
	    for (int i = 0; i < _m * _n; ++i) {
		_data[i] = C._data[i];
	    }
	}
    }
};

template <class F> inline std::ostream& operator<<(std::ostream& os,
                                                   const NumMat<F>& mat) {
#ifndef RELEASE
    CallStackEntry entry("operator<<");
#endif
    os << mat.m() << " " << mat.n() << std::endl;
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for (int i = 0; i < mat.m(); ++i) {
        for (int j = 0; j < mat.n(); ++j)
            os << " " << mat(i,j);
        os << std::endl;
    }
    return os;
}
template <class F> inline void setvalue(NumMat<F>& M, F val) {
#ifndef RELEASE
    CallStackEntry entry("setvalue");
#endif
    for (int i = 0; i < M.m(); ++i) {
        for (int j = 0; j < M.n(); ++j) {
            M(i,j) = val;
        }
    }
}
template <class F> inline double energy(NumMat<F>& M) {
#ifndef RELEASE
    CallStackEntry entry("energy");
#endif
    double sum = 0;
    for (int i = 0; i < M.m(); ++i) {
        for (int j = 0; j < M.n(); ++j)  {
            sum += abs(M(i,j) * M(i,j));
        }
    }
    return sum;
}

typedef NumMat<bool>   BolNumMat;
typedef NumMat<int>    IntNumMat;
typedef NumMat<double> DblNumMat;
typedef NumMat<cpx>    CpxNumMat;

#endif
