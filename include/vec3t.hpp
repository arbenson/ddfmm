#ifndef  _VEC3T_HPP_
#define  _VEC3T_HPP_

#include "commoninc.hpp"

using std::istream;
using std::ostream;
using std::min;
using std::max;
using std::abs;

// Common VECtor Template
template <class F>
class Vec3T {
private:
    F _v[3];
public:

    enum { X=0, Y=1, Z=2 };

    //------------CONSTRUCTOR AND DESTRUCTOR
    Vec3T() {
        _v[0] = F(0);
        _v[1] = F(0);
        _v[2] = F(0);
    }

    Vec3T(const F* f) {
        _v[0] = f[0];
        _v[1] = f[1];
        _v[2] = f[2];
    }

    Vec3T(F a, F b, F c)   {
        _v[0] = a;
        _v[1] = b;
        _v[2] = c;
    }

    Vec3T(const Vec3T& c) {
        _v[0] = c._v[0];
        _v[1] = c._v[1];
        _v[2] = c._v[2];
    }

    ~Vec3T() {}

    //------------POINTER and ACCESS
    operator F*()             { return &_v[0]; }
    operator const F*() const { return &_v[0]; }
    F* array()                { return &_v[0]; }  //access array
    F& operator()(int i)             { assert(i<3); return _v[i]; }
    const F& operator()(int i) const { assert(i<3); return _v[i]; }
    F& operator[](int i)             { assert(i<3); return _v[i]; }
    const F& operator[](int i) const { assert(i<3); return _v[i]; }
    //------------ASSIGN
    Vec3T& operator= ( const Vec3T& c ) { _v[0] =c._v[0]; _v[1] =c._v[1]; _v[2] =c._v[2]; return *this; }
    Vec3T& operator+=( const Vec3T& c ) { _v[0]+=c._v[0]; _v[1]+=c._v[1]; _v[2]+=c._v[2]; return *this; }
    Vec3T& operator-=( const Vec3T& c ) { _v[0]-=c._v[0]; _v[1]-=c._v[1]; _v[2]-=c._v[2]; return *this; }
    Vec3T& operator*=( const F& s )     { _v[0]*=s;       _v[1]*=s;       _v[2]*=s;       return *this; }
    Vec3T& operator/=( const F& s )     { _v[0]/=s;       _v[1]/=s;       _v[2]/=s;       return *this; }
    //-----------LENGTH...
    F l1( void )     const  { F sum=F(0); for(int i=0; i<3; i++) sum=sum+abs(_v[i]); return sum; }
    F linfty( void ) const  { F cur=F(0); for(int i=0; i<3; i++) cur=max(cur,abs(_v[i])); return cur; }
    F l2( void )     const  { F sum=F(0); for(int i=0; i<3; i++) sum=sum+_v[i]*_v[i]; return sqrt(sum); }
};

//-----------LEX COMPARE
template <class F> inline bool operator==(const Vec3T<F>& a, const Vec3T<F>& b) {
    return (a[0]==b[0] && a[1]==b[1] && a[2]==b[2]);
}
template <class F> inline bool operator!=(const Vec3T<F>& a, const Vec3T<F>& b) {
    return !(a==b);
}
template <class F> inline bool operator> (const Vec3T<F>& a, const Vec3T<F>& b) {
    for(int i=0; i<3; i++) {
        if (a[i] > b[i]) {
            return true;
        } else if (a[i] < b[i]) {
            return false;
        }
    }
    return false;
}
template <class F> inline bool operator< (const Vec3T<F>& a, const Vec3T<F>& b) {
    for(int i = 0; i < 3; i++) {
        if (a[i] < b[i]) {
            return true;
        } else if (a[i] > b[i]) {
            return false;
        }
    }
    return false;
}
template <class F> inline bool operator>=(const Vec3T<F>& a, const Vec3T<F>& b) {
    for(int i = 0; i < 3; i++) {
        if (a[i] > b[i]) {
            return true;
        }       else if (a[i] < b[i]) {
            return false;
        }
    }
    return true;
}
template <class F> inline bool operator<=(const Vec3T<F>& a, const Vec3T<F>& b) {
    for(int i = 0; i < 3; i++) {
        if (a[i] < b[i]) {
            return true;
        } else if (a[i] > b[i]) {
            return false;
        }
    }
    return true;
}

//-----------NUMERICAL OPS
template <class F> inline Vec3T<F> operator- (const Vec3T<F>& a) {
    return Vec3T<F>(-a[0], -a[1], -a[2]);
}
template <class F> inline Vec3T<F> operator+ (const Vec3T<F>& a, const Vec3T<F>& b) {
    return Vec3T<F>(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
template <class F> inline Vec3T<F> operator- (const Vec3T<F>& a, const Vec3T<F>& b) {
    return Vec3T<F>(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
template <class F> inline Vec3T<F> operator* (F scl, const Vec3T<F>& a) {
    return Vec3T<F>(scl * a[0], scl * a[1], scl * a[2]);
}
template <class F> inline Vec3T<F> operator* (const Vec3T<F>& a, F scl) {
    return Vec3T<F>(scl * a[0], scl * a[1], scl * a[2]);
}
template <class F> inline Vec3T<F> operator/ (const Vec3T<F>& a, F scl) {
    return Vec3T<F>(a[0] / scl, a[1] / scl, a[2] / scl);
}
template <class F> inline F operator* (const Vec3T<F>& a, const Vec3T<F>& b) {
    F sum = F(0); 
    sum += a(0) * b(0);
    sum += a(1) * b(1);
    sum += a(2) * b(2);
    return sum;
}
template <class F> inline F dot (const Vec3T<F>& a, const Vec3T<F>& b) {
    return a*b;
}
template <class F> inline Vec3T<F> operator^ (const Vec3T<F>& a, const Vec3T<F>& b) {
    return Vec3T<F>(a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0));
}
template <class F> inline Vec3T<F> cross (const Vec3T<F>& a, const Vec3T<F>& b) { 
    return a^b; 
}

//-------------ew NUMERICAL OPS
template <class F> inline Vec3T<F> ewmin(const Vec3T<F>& a, const Vec3T<F>& b) {
    return Vec3T<F>(min(a[0], b[0]), min(a[1], b[1]), min(a[2], b[2]));
}
template <class F> inline Vec3T<F> ewmax(const Vec3T<F>& a, const Vec3T<F>& b) {
    return Vec3T<F>(max(a[0], b[0]), max(a[1], b[1]), max(a[2], b[2]));
}
template <class F> inline Vec3T<F> ewabs(const Vec3T<F>& a) {
    return Vec3T<F>(abs(a[0]), abs(a[1]), abs(a[2]));
}
template <class F> inline Vec3T<F> ewmul(const Vec3T<F>&a, const Vec3T<F>& b) {
    return Vec3T<F>(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
template <class F> inline Vec3T<F> ewdiv(const Vec3T<F>&a, const Vec3T<F>& b) { 
    return Vec3T<F>(a[0] / b[0], a[1] / b[1], a[2] / b[2]);
}
template <class F> inline Vec3T<F> ewrnd(const Vec3T<F>&a) { //round
    return Vec3T<F>(round(a[0]), round(a[1]), round(a[2]));
}

//---------------INOUT
template <class F> istream& operator>>(istream& is, Vec3T<F>& a) {
    for(int i=0; i<3; i++)
        is >> a[i];
    return is;
}
template <class F> ostream& operator<<(ostream& os, const Vec3T<F>& a) { 
    for(int i=0; i<3; i++)
        os << a[i] << " ";
    return os;
}

//---------------------------------------------------------
/// MOST COMMONLY USED
typedef Vec3T<double> Point3;
typedef Vec3T<int>    Index3;

#endif

