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
#ifndef _SERIALIZE_HPP_
#define _SERIALIZE_HPP_

#include "commoninc.hpp"
#include "wave3d.hpp"

template<class T, class S> int serialize(const std::pair<T,S>& val,
                                         std::ostream& os,
					 const std::vector<int>& mask);
template<class T, class S> int deserialize(std::pair<T,S>& val,
                                           std::istream& is,
					   const std::vector<int>& mask);

//-------------------
//char
inline int serialize(const char& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    os.write((char*)&val, sizeof(char));
    return 0;
}

inline int deserialize(char& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    is.read((char*)&val, sizeof(char));
    return 0;
}

//-------------------
//int
inline int serialize(const int& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    os.write((char*)&val, sizeof(int));
    return 0;
}

inline int deserialize(int& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    is.read((char*)&val, sizeof(int));
    return 0;
}

//-------------------
//double
inline int serialize(const double& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    os.write((char*)&val, sizeof(double));
    return 0;
}

inline int deserialize(double& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    is.read((char*)&val, sizeof(double));
    return 0;
}

//-------------------
//cpx
inline int serialize(const cpx& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    os.write((char*)&val, sizeof(cpx));
    return 0;
}

inline int deserialize(cpx& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    is.read((char*)&val, sizeof(cpx));
    return 0;
}


//-------------------
//Index3
inline int serialize(const Index3& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    os.write((char*)&(val[0]), 3*sizeof(int));
    return 0;
}

inline int deserialize(Index3& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    is.read((char*)&(val[0]), 3*sizeof(int));
    return 0;
}

//-------------------
//Point3
inline int serialize(const Point3& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    os.write((char*)&(val[0]), 3*sizeof(double));
    return 0;
}

inline int deserialize(Point3& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    is.read((char*)&(val[0]), 3*sizeof(double));
    return 0;
}

//-------------------
//std::vector
template<class T>
int serialize(const std::vector<T>& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int sz = val.size();
    os.write((char*)&sz, sizeof(int));
    for (int k = 0; k < sz; k++)
        serialize(val[k], os, mask);
    return 0;
}

template<class T>
int deserialize(std::vector<T>& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int sz;
    is.read((char*)&sz, sizeof(int));
    val.resize(sz);
    for (int k = 0; k < sz; k++)
        deserialize(val[k], is, mask);
    return 0;
}

//-------------------
//set
template<class T>
int serialize(const std::set<T>& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int sz = val.size();
    os.write((char*)&sz, sizeof(int));
    for (typename std::set<T>::const_iterator mi = val.begin(); mi != val.end(); mi++)
        serialize((*mi), os, mask);
    return 0;
}

template<class T>
int deserialize(std::set<T>& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    val.clear();
    int sz;
    is.read((char*)&sz, sizeof(int));
    for(int k=0; k<sz; k++) {
        T t; deserialize(t, is, mask);
        val.insert(t);
    }
    return 0;
}

//-------------------
//map
template<class T, class S>
int serialize(const std::map<T, S>& val, std::ostream& os,
              const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int sz = val.size();
    os.write((char*)&sz, sizeof(int));
    for (typename std::map<T, S>::const_iterator mi = val.begin(); mi != val.end();
         mi++) {
        serialize(mi->first, os, mask);
        serialize(mi->second, os, mask);
    }
    return 0;
}

template<class T, class S>
int deserialize(std::map<T, S>& val, std::istream& is,
                const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    val.clear();
    int sz;
    is.read((char*)&sz, sizeof(int));
    for (int k = 0; k < sz; k++) {
        T t;
	deserialize(t, is, mask);
        S s;
	deserialize(s, is, mask);
        val[t] = s;
    }
    return 0;
}

//-------------------
//pair
template<class T, class S>
int serialize(const std::pair<T,S>& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    serialize(val.first, os, mask);
    serialize(val.second, os, mask);
    return 0;
}

template<class T, class S>
int deserialize(std::pair<T,S>& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    deserialize(val.first, is, mask);
    deserialize(val.second, is, mask);
    return 0;
}


//-------------------
//BolNumVec
inline int serialize(const BolNumVec& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    os.write((char*)&m, sizeof(int));
    os.write((char*)(val.data()), m*sizeof(bool));
    return 0;
}

inline int deserialize(BolNumVec& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    is.read((char*)&m, sizeof(int));
    val.resize(m);
    is.read((char*)(val.data()), m*sizeof(bool));
    return 0;
}

//-------------------
//BolNumMat
inline int serialize(const BolNumMat& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    int n = val.n();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    os.write((char*)(val.data()), m*n*sizeof(bool));
    return 0;
}

inline int deserialize(BolNumMat& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    int n;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    val.resize(m,n);
    is.read((char*)(val.data()), m*n*sizeof(bool));
    return 0;
}

//-------------------
//BolNumTns
inline int serialize(const BolNumTns& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();  int n = val.n();  int p = val.p();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    os.write((char*)&p, sizeof(int));
    os.write((char*)(val.data()), m*n*p*sizeof(bool));
    return 0;
}

inline int deserialize(BolNumTns& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m,n,p;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    is.read((char*)&p, sizeof(int));
    val.resize(m,n,p);
    is.read((char*)(val.data()), m*n*p*sizeof(bool));
    return 0;
}


//-------------------
//IntNumVec
inline int serialize(const IntNumVec& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    os.write((char*)&m, sizeof(int));
    os.write((char*)(val.data()), m*sizeof(int));
    return 0;
}

inline int deserialize(IntNumVec& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    is.read((char*)&m, sizeof(int));
    val.resize(m);
    is.read((char*)(val.data()), m*sizeof(int));
    return 0;
}

//-------------------
//IntNumMat
inline int serialize(const IntNumMat& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    int n = val.n();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    os.write((char*)(val.data()), m*n*sizeof(int));
    return 0;
}

inline int deserialize(IntNumMat& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    int n;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    val.resize(m,n);
    is.read((char*)(val.data()), m*n*sizeof(int));
    return 0;
}

//-------------------
//IntNumTns
inline int serialize(const IntNumTns& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();  int n = val.n();  int p = val.p();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    os.write((char*)&p, sizeof(int));
    os.write((char*)(val.data()), m*n*p*sizeof(int));
    return 0;
}

inline int deserialize(IntNumTns& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m,n,p;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    is.read((char*)&p, sizeof(int));
    val.resize(m,n,p);
    is.read((char*)(val.data()), m*n*p*sizeof(int));
    return 0;
}


//-------------------
//DblNumVec
inline int serialize(const DblNumVec& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    os.write((char*)&m, sizeof(int));
    os.write((char*)(val.data()), m*sizeof(double));
    return 0;
}

inline int deserialize(DblNumVec& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    is.read((char*)&m, sizeof(int));
    val.resize(m);
    is.read((char*)(val.data()), m*sizeof(double));
    return 0;
}

//-------------------
//DblNumMat
inline int serialize(const DblNumMat& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    int n = val.n();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    os.write((char*)(val.data()), m*n*sizeof(double));
    return 0;
}

inline int deserialize(DblNumMat& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    int n;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    val.resize(m,n);
    is.read((char*)(val.data()), m*n*sizeof(double));
    return 0;
}

//-------------------
//DblNumTns
inline int serialize(const DblNumTns& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();  int n = val.n();  int p = val.p();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    os.write((char*)&p, sizeof(int));
    os.write((char*)(val.data()), m*n*p*sizeof(double));
    return 0;
}

inline int deserialize(DblNumTns& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m,n,p;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    is.read((char*)&p, sizeof(int));
    val.resize(m,n,p);
    is.read((char*)(val.data()), m*n*p*sizeof(double));
    return 0;
}

//-------------------
//CpxNumVec
inline int serialize(const CpxNumVec& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    os.write((char*)&m, sizeof(int));
    os.write((char*)(val.data()), m*sizeof(cpx));
    return 0;
}

inline int deserialize(CpxNumVec& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    is.read((char*)&m, sizeof(int));
    val.resize(m);
    is.read((char*)(val.data()), m*sizeof(cpx));
    return 0;
}

//-------------------
//CpxNumMat
inline int serialize(const CpxNumMat& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    int n = val.n();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    os.write((char*)(val.data()), m*n*sizeof(cpx));
    return 0;
}

inline int deserialize(CpxNumMat& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    int n;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    val.resize(m,n);
    is.read((char*)(val.data()), m*n*sizeof(cpx));
    return 0;
}

//-------------------
//CpxNumTns
inline int serialize(const CpxNumTns& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();  int n = val.n();  int p = val.p();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    os.write((char*)&p, sizeof(int));
    os.write((char*)(val.data()), m*n*p*sizeof(cpx));
    return 0;
}

inline int deserialize(CpxNumTns& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m,n,p;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    is.read((char*)&p, sizeof(int));
    val.resize(m,n,p);
    is.read((char*)(val.data()), m*n*p*sizeof(cpx));
    return 0;
}

//-------------------
//NumVec
template<class T>
int inline serialize(const NumVec<T>& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    os.write((char*)&m, sizeof(int));
    for(int i=0; i<m; i++)
        serialize(val(i), os, mask);
    return 0;
}
template<class T>
int inline deserialize(NumVec<T>& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    is.read((char*)&m, sizeof(int));
    val.resize(m);
    for (int i = 0; i < m; i++)
        deserialize(val(i), is, mask);
    return 0;
}

//-------------------
//NumMat
template<class T>
int inline serialize(const NumMat<T>& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    int n = val.n();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            serialize(val(i,j), os, mask);
    return 0;
}
template<class T>
int inline deserialize(NumMat<T>& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    int n;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    val.resize(m,n);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            deserialize(val(i,j), is, mask);
    return 0;
}

//-------------------
//NumTns
template<class T>
int inline serialize(const NumTns<T>& val, std::ostream& os, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    int m = val.m();
    int n = val.n();
    int p = val.p();
    os.write((char*)&m, sizeof(int));
    os.write((char*)&n, sizeof(int));
    os.write((char*)&p, sizeof(int));
    for (int k = 0; k < p; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                serialize(val(i,j,k), os, mask);
    return 0;
}
template<class T>
int inline deserialize(NumTns<T>& val, std::istream& is, const std::vector<int>& mask)
{
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    int m;
    int n;
    int p;
    is.read((char*)&m, sizeof(int));
    is.read((char*)&n, sizeof(int));
    is.read((char*)&p, sizeof(int));
    val.resize(m,n,p);
    for (int k = 0; k < p; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                deserialize(val(i,j,k), is, mask);
    return 0;
}


#endif
