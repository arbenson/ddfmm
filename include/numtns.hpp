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
#ifndef _NUMTNS_HPP_
#define _NUMTNS_HPP_

#include "nummat.hpp"

template <class F>
class NumTns {
public:
    NumTns(int m=0, int n=0, int p=0): _m(m), _n(n), _p(p), _owndata(true) {
#ifndef RELEASE
        CallStackEntry entry("NumTns::NumTns");
#endif
        allocate();
    }

    NumTns(int m, int n, int p, bool owndata, F* data): _m(m), _n(n), _p(p),
                                                        _owndata(owndata) {
#ifndef RELEASE
        CallStackEntry entry("NumTns::NumTns");
#endif
        if (_owndata) {
            allocate();
            if (ValidDimensions()) {
                for (int i = 0; i < _m * _n * _p; ++i) {
                    _data[i] = data[i];
                }
            }
        } else {
            _data = data;
        }
    }

    NumTns(const NumTns& C): _m(C._m), _n(C._n), _p(C._p),
                             _owndata(C._owndata) {
#ifndef RELEASE
        CallStackEntry entry("NumTns::NumTns");
#endif
        if (_owndata) {
            allocate();
            fill(C);
        } else {
            _data = C._data;
        }
    }

    ~NumTns() { 
#ifndef RELEASE
        CallStackEntry entry("NumTns::~NumTns");
#endif
        if (_owndata)
            deallocate();
    }

    NumTns& operator=(const NumTns& C) {
#ifndef RELEASE
        CallStackEntry entry("NumTns::operator=");
#endif
        if (_owndata)
            deallocate();
        _m = C._m;
        _n = C._n;
        _p = C._p;
        _owndata = C._owndata;
        if (_owndata) {
            allocate();
            fill(C);
        } else {
            _data = C._data;
        }
        return *this;
    }

    void resize(int m, int n, int p)  {
#ifndef RELEASE
        CallStackEntry entry("NumTns::resize");
#endif
        assert( _owndata );
        if (_m != m || _n != n || _p != p) {
            deallocate();
            _m = m;
            _n = n;
            _p = p;
            allocate();
        }
    }

    const F& operator()(int i, int j, int k) const  {
#ifndef RELEASE
        CallStackEntry entry("NumTns::operator()");
#endif
        assert( i >= 0 && i < _m && j >= 0 && j < _n && k >= 0 && k < _p);
        return _data[i + j * _m + k * _m * _n];
    }
    F& operator()(int i, int j, int k)  {
#ifndef RELEASE
        CallStackEntry entry("NumTns::operator()");
#endif
        assert( i >= 0 && i < _m && j >= 0 && j < _n && k >= 0 && k < _p);
        return _data[i + j * _m + k * _m * _n];
    }
    
    F* data() const { return _data; }
    int m() const { return _m; }
    int n() const { return _n; }
    int p() const { return _p; }

private:
    int _m, _n, _p;
    bool _owndata;
    F* _data;

    inline bool ValidDimensions() const { return _m > 0 && _n > 0 && _p > 0; }
    
    void allocate() {
#ifndef RELEASE
        CallStackEntry entry("NumTns::allocate");
#endif
        if (ValidDimensions()) {
            _data = new F[_m * _n * _p];
            assert( _data != NULL );
        } else {
            _data = NULL;
        }
    }

    void deallocate() {
#ifndef RELEASE
        CallStackEntry entry("NumTns::deallocate");
#endif
        if (ValidDimensions()) {
            delete[] _data;
            _data = NULL;
        }
    }

    void fill(const NumTns& C) {
#ifndef RELEASE
        CallStackEntry entry("NumTns::fill");
#endif
        if (ValidDimensions()) {
            for (int i = 0; i < _m * _n * _p; ++i) {
                _data[i] = C._data[i];
            }
        }
    }
};

template <class F> inline std::ostream& operator<<(std::ostream& os,
                                                   const NumTns<F>& tns) {
#ifndef RELEASE
    CallStackEntry entry("operator<<");
#endif
    os << tns.m() << " " << tns.n() << " " << tns.p() << std::endl;
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for (int i = 0; i < tns.m(); ++i) {
        for (int j = 0; j < tns.n(); ++j) {
            for (int k = 0; k < tns.p(); ++k) {
                os << " " << tns(i,j,k);
            }
            os << std::endl;
        }
        os << std::endl;
    }
    return os;
}

template <class F> inline void setvalue(NumTns<F>& T, F val) {
#ifndef RELEASE
    CallStackEntry entry("setvalue");
#endif
    for (int i = 0; i < T.m(); ++i) {
        for (int j = 0; j < T.n(); ++j) {
            for (int k = 0; k < T.p(); ++k) {
                T(i,j,k) = val;
            }
        }
    }
  return;
}

template <class F> inline double energy(NumTns<F>& T) {
#ifndef RELEASE
    CallStackEntry entry("energy");
#endif
  double sum = 0;
  for (int i = 0; i < T.m(); ++i) {
      for (int j = 0; j < T.n(); ++j) {
          for (int k = 0; k < T.p(); ++k) {
              sum += abs(T(i,j,k) * T(i,j,k));
          }
      }
  }
  return sum;
}

template <class F> inline double NumTnsSum(NumTns<F>& T) {
#ifndef RELEASE
    CallStackEntry entry("energy");
#endif
  double sum = 0;
  for (int i = 0; i < T.m(); ++i) {
      for (int j = 0; j < T.n(); ++j) {
          for (int k = 0; k < T.p(); ++k) {
              sum += T(i,j,k);
          }
      }
  }
  return sum;
}

typedef NumTns<bool>   BolNumTns;
typedef NumTns<int>    IntNumTns;
typedef NumTns<double> DblNumTns;
typedef NumTns<cpx>    CpxNumTns;

#endif
