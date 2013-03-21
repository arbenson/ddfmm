/* Kernel Independent Fast Multipole Method
   Copyright (C) 2004 Lexing Ying, New York University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */
#ifndef  _VEC2T_HPP_
#define  _VEC2T_HPP_

#include "commoninc.hpp"

using std::istream;
using std::ostream;
using std::min;
using std::max;
using std::abs;

//-----------------------------------------------------------------------------------------------------
///Common VECtor Template
template <class F>
class Vec2T {
private:
  F _v[2];
public:
  enum{ X=0, Y=1 };
  //------------CONSTRUCTOR AND DESTRUCTOR 
  Vec2T()              { _v[0]=F(0);    _v[1]=F(0); }
  Vec2T(F f)           { _v[0]=f;       _v[1]=f; }
  Vec2T(const F* f)    { _v[0]=f[0];    _v[1]=f[1]; }
  Vec2T(F a,F b)       { _v[0]=a;       _v[1]=b; }
  Vec2T(const Vec2T& c){ _v[0]=c._v[0]; _v[1]=c._v[1]; }
  ~Vec2T() {}
  //------------POINTER and ACCESS
  operator F*()             { return &_v[0]; }
  operator const F*() const { return &_v[0]; }
  F* array()                { return &_v[0]; }  //access array
  F& operator()(int i)             { assert(i<2); return _v[i]; }
  const F& operator()(int i) const { assert(i<2); return _v[i]; }
  F& operator[](int i)             { assert(i<2); return _v[i]; }
  const F& operator[](int i) const { assert(i<2); return _v[i]; }
  F& x()             { return _v[0];}
  F& y()             { return _v[1];}
  const F& x() const { return _v[0];}
  const F& y() const { return _v[1];}
  //------------ASSIGN
  Vec2T& operator= ( const Vec2T& c ) { _v[0] =c._v[0]; _v[1] =c._v[1]; return *this; }
  Vec2T& operator+=( const Vec2T& c ) { _v[0]+=c._v[0]; _v[1]+=c._v[1]; return *this; }
  Vec2T& operator-=( const Vec2T& c ) { _v[0]-=c._v[0]; _v[1]-=c._v[1]; return *this; }
  Vec2T& operator*=( const F& s )     { _v[0]*=s;       _v[1]*=s;       return *this; }
  Vec2T& operator/=( const F& s )     { _v[0]/=s;       _v[1]/=s;       return *this; }
  //-----------LENGTH...
  F l1( void )     const  { F sum=F(0); for(int i=0; i<2; i++) sum=sum+abs(_v[i]); return sum; }
  F linfty( void ) const  { F cur=F(0); for(int i=0; i<2; i++) cur=max(cur,abs(_v[i])); return cur; }
  F l2( void )     const  { F sum=F(0); for(int i=0; i<2; i++) sum=sum+_v[i]*_v[i]; return sqrt(sum); }
  F length( void ) const  { return l2(); }
  Vec2T dir( void )    const  { F a=l2(); return (*this)/a; }
};

//-----------BOOLEAN OPS
template <class F> inline bool operator==(const Vec2T<F>& a, const Vec2T<F>& b) {
  bool res = true;  for(int i=0; i<2; i++)   res = res && (a(i)==b(i));  return res;
}
template <class F> inline bool operator!=(const Vec2T<F>& a, const Vec2T<F>& b) {
  return !(a==b);
}
template <class F> inline bool operator> (const Vec2T<F>& a, const Vec2T<F>& b) {
  bool res = true;  for(int i=0; i<2; i++)   res = res && (a(i)> b(i));  return res; 
}
template <class F> inline bool operator< (const Vec2T<F>& a, const Vec2T<F>& b) {
  bool res = true;  for(int i=0; i<2; i++)   res = res && (a(i)< b(i));  return res; 
}
template <class F> inline bool operator>=(const Vec2T<F>& a, const Vec2T<F>& b) {
  bool res = true;  for(int i=0; i<2; i++)	res = res && (a(i)>=b(i));  return res; 
}
template <class F> inline bool operator<=(const Vec2T<F>& a, const Vec2T<F>& b) {
  bool res = true;  for(int i=0; i<2; i++)   res = res && (a(i)<=b(i));  return res; 
}

//-----------NUMERICAL OPS
template <class F> inline Vec2T<F> operator- (const Vec2T<F>& a) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = -a[i]; return r;
}
template <class F> inline Vec2T<F> operator+ (const Vec2T<F>& a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]+b[i]; return r; 
}
template <class F> inline Vec2T<F> operator- (const Vec2T<F>& a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]-b[i]; return r;
}
template <class F> inline Vec2T<F> operator* (F scl, const Vec2T<F>& a) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec2T<F> operator* (const Vec2T<F>& a, F scl) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec2T<F> operator/ (const Vec2T<F>& a, F scl) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]/scl;  return r;
}

template <class F> inline F operator* (const Vec2T<F>& a, const Vec2T<F>& b) {
  F sum=F(0); for(int i=0; i<2; i++) sum=sum+a(i)*b(i); return sum;
}
template <class F> inline F dot       (const Vec2T<F>& a, const Vec2T<F>& b) {
  return a*b;
}
//-------------ew OPS
template <class F> inline Vec2T<F> min(const Vec2T<F>& a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = min(a[i], b[i]); return r;
}
template <class F> inline Vec2T<F> max(const Vec2T<F>& a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = max(a[i], b[i]); return r;
}
template <class F> inline Vec2T<F> abs(const Vec2T<F>& a) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = abs(a[i]); return r;
}
template <class F> inline Vec2T<F> ewmul(const Vec2T<F>&a, const Vec2T<F>& b) {
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]*b[i]; return r;
}
template <class F> inline Vec2T<F> ewdiv(const Vec2T<F>&a, const Vec2T<F>& b) { 
  Vec2T<F> r;  for(int i=0; i<2; i++) r[i] = a[i]/b[i]; return r;
}
//---------------INOUT
template <class F> istream& operator>>(istream& is, Vec2T<F>& a) {
  for(int i=0; i<2; i++) is>>a[i]; return is;
}
template <class F> ostream& operator<<(ostream& os, const Vec2T<F>& a) { 
  for(int i=0; i<2; i++) os<<a[i]<<" "; return os;
}

//---------------------------------------------------------
/// MOST COMMONLY USED
typedef Vec2T<double> Point2;
typedef Vec2T<int>    Index2;



#endif
