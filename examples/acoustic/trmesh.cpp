#include "trmesh.hpp"

using std::map;
using std::ifstream;
using std::istringstream;

// ----------------------------------------------------------------------
int TrMesh::setup(vector<Point3>& points, vector<Index3>& coords)
{
  //nv nf ne
  int nv = points.size();
  int nf = coords.size();
  int ne = coords.size()*3;
  _edges.resize(ne);  _verts.resize(nv);  _faces.resize(nf);
  //_faces, the easiest
  for(int F=0; F<nf; F++) {
    _faces[F].SE = 3*F;
  }
  //_edges,
  map<intpair,int> VV2E;
  int ec=0; //edge count
  for(int F=0; F<nf; F++) {
    for(int e=0; e<3; e++) {
      int fm = coords[F][fromvno(e)];
      int to = coords[F][gotovno(e)];
      _edges[ec+e].SV = fm;
      _edges[ec+e].CE = -1; //tmp
      _edges[ec+e].F = F;
      _edges[ec+e].PE = ec+preveno(e);
      _edges[ec+e].NE = ec+nexteno(e);
      VV2E[intpair(fm,to)] = ec+e;
    }
    ec += 3;
  }
  CHECK_TRUE(ec==ne);
  ec = 0;
  for(int F=0; F<nf; F++) {
    for(int e=0; e<3; e++) {
      int fm = coords[F][fromvno(e)];
      int to = coords[F][gotovno(e)];
      map<intpair,int>::iterator mi = VV2E.find(intpair(to,fm));
      if(mi!=VV2E.end())
	_edges[ec+e].CE = (*mi).second;
    }
    ec += 3;
  }
  CHECK_TRUE(ec==ne);
  //_verts,
  for(int V=0; V<nv; V++)
    _verts[V].SE = -1; //set invalid
  for(int E=0; E<ne; E++) {
    int V = _edges[E].SV; //start v
    if(_verts[V].SE==-1) { //do only once
      bool start = true;
      int ME = E;
      while( _edges[ME].CE !=-1 && (start==true||ME!=E)  ) {
	start = false;
	ME = _edges[_edges[ME].CE].NE;
      }
      _verts[V].SE = ME;
    }
  }
  //points
  _points = points;
  return (0);
}

// ----------------------------------------------------------------------
int TrMesh::compute_dihedral(vector<double>& dihs)
{
  int ne = _edges.size();
  dihs.resize(ne);
  for(int e=0; e<ne; e++) {
    int e0 = e;
    int e1 = _edges[e0].CE;
    int v0 = _edges[e0].SV;
    int v1 = _edges[e1].SV;
    int e2 = _edges[e0].PE;
    int e3 = _edges[e1].PE;
    int v2 = _edges[e2].SV;
    int v3 = _edges[e3].SV;
    Point3 p0 = _points[v0];
    Point3 p1 = _points[v1];
    Point3 p2 = _points[v2];
    Point3 p3 = _points[v3];
    Point3 n0 = cross(p0-p2, p1-p0);    n0 = n0/n0.l2();
    Point3 n1 = cross(p1-p3, p0-p1);    n1 = n1/n1.l2();
    Point3 t0 = cross(n0, p1-p0);    t0 = t0/t0.l2();
    Point3 t1 = cross(n1, p0-p1);    t1 = t1/t1.l2();
    Point3 zd = p0-p1;    zd = zd/zd.l2();
    double cosang = dot(t0,t1);
    double sinang = dot(zd, cross(t0,t1));
    double tmp = atan2(sinang, cosang);
    if(tmp<0) tmp = tmp+2*M_PI;
    dihs[e0] = tmp; //angles
  }
  return 0;
}

// ----------------------------------------------------------------------
int TrMesh::compute_interior(vector<double>& itrs)
{
  vector<double> dihs;
  SAFE_FUNC_EVAL( compute_dihedral(dihs) );
  //
  int nv = _verts.size();
  itrs.resize(nv);
  for(int v=0; v<nv; v++) {
    //1. get one ring of edges
    vector<int> es;
    int estt = _verts[v].SE;
    es.push_back(estt);
    int enow = _edges[_edges[estt].PE].CE;
    while(enow!=estt) {
      es.push_back(enow);
      enow = _edges[_edges[enow].PE].CE;
    }
    double sum = 0;
    for(int k=0; k<es.size(); k++) {
      sum = sum + dihs[ es[k] ];
    }
    sum = sum - M_PI*(es.size()-2);
    itrs[v] = sum;
  }
  return 0;
}

// ----------------------------------------------------------------------
int TrMesh::compute_area(vector<double>& ares)
{
  int nf = _faces.size();
  ares.resize(nf);
  for(int f=0; f<nf; f++) {
    vector<int> vs;
    F2Vs(f, vs);
    Point3 v0 = _points[vs[0]];
    Point3 v1 = _points[vs[1]];
    Point3 v2 = _points[vs[2]];
    ares[f] = (cross(v1-v0,v2-v0).l2())/2.0;
  }
  return 0;
}

