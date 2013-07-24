#ifndef _NUMVEC_HPP_
#define _NUMVEC_HPP_

#include "commoninc.hpp"

using std::istream;
using std::ostream;
using std::ios_base;
using std::endl;

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

template <class F> inline ostream& operator<<( ostream& os, const NumVec<F>& vec)
{
#ifndef RELEASE
    CallStackEntry entry("operator<<");
#endif
    os << vec.m() << endl;
    os.setf(ios_base::scientific, ios_base::floatfield);
    for(int i = 0; i < vec.m(); i++) {
        os << " " << vec(i);
    }
    os << endl;
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


