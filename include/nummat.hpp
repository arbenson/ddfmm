#ifndef _NUMMAT_HPP_
#define _NUMMAT_HPP_

#include "numvec.hpp"

template <class F>
class NumMat
{
public:
    NumMat(int m=0, int n=0): _m(m), _n(n), _owndata(true) {
#ifndef RELEASE
        CallStackEntry entry("NumMat::NumMat");
#endif
	allocate();
    }

    NumMat(int m, int n, bool owndata, F* data): _m(m), _n(n), _owndata(owndata) {
#ifndef RELEASE
        CallStackEntry entry("NumMat::NumMat");
#endif
        if (_owndata) {
	    allocate();
            if (checkDimensions()) {
                for(int i = 0; i < _m * _n; i++) {
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
            if (checkDimensions()) {
                delete[] _data;
                _data = NULL;
            }
            _m = m;
            _n = n;
	    allocate();
        }
    }

    const F& operator()(int i, int j) const  { 
#ifndef RELEASE
        CallStackEntry entry("NumMat::operator()");
#endif
        assert( i >= 0 && i < _m && j >= 0 && j <_ n );
        return _data[i + j * _m];
    }

    F& operator()(int i, int j)  { 
#ifndef RELEASE
        CallStackEntry entry("NumMat::operator()");
#endif
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

    bool checkDimensions() const { return _m > 0 && _n > 0; }

    void allocate() {
#ifndef RELEASE
        CallStackEntry entry("NumMat::allocate");
#endif
	if (checkDimensions()) {
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
	if (checkDimensions()) {
	    delete[] _data;
	    _data = NULL;
	}
    }

    void fill(const NumMat& C) {
#ifndef RELEASE
        CallStackEntry entry("NumMat::fill");
#endif
	if (checkDimensions()) {
	    for (int i = 0; i < _m * _n; i++) {
		_data[i] = C._data[i];
	    }
	}
    }
};

template <class F> inline ostream& operator<<( ostream& os, const NumMat<F>& mat)
{
#ifndef RELEASE
    CallStackEntry entry("operator<<");
#endif
    os << mat.m() << " " << mat.n() << endl;
    os.setf(ios_base::scientific, ios_base::floatfield);
    for(int i = 0; i < mat.m(); i++) {
        for(int j = 0; j < mat.n(); j++)
            os << " " << mat(i,j);
        os << endl;
    }
    return os;
}
template <class F> inline void setvalue(NumMat<F>& M, F val)
{
#ifndef RELEASE
    CallStackEntry entry("setvalue");
#endif
    for (int i = 0; i < M.m(); i++) {
        for (int j = 0; j < M.n(); j++) {
            M(i,j) = val;
        }
    }
}
template <class F> inline double energy(NumMat<F>& M)
{
#ifndef RELEASE
    CallStackEntry entry("energy");
#endif
    double sum = 0;
    for (int i = 0; i < M.m(); i++) {
        for (int j = 0; j < M.n(); j++)  {
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
