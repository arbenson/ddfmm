#ifndef _TRMESH_HPP_
#define _TRMESH_HPP_

#include "comobject.hpp"
#include "vec3t.hpp"

#include <string>

using std::vector;
using std::pair;

//----------------------
//triangle mesh
//----------------------
class TrMesh: public ComObject
{
public:
  typedef pair<int,int> intpair;
  //-----------
  //Topology related data types
  struct Edge {
    int SV; //starting vert
    int CE; //conjugate(neighbor) edge
    int F;  //face
    int PE; //prev edge
    int NE; //next edge
  };
  struct Vert {
    int SE; //starting edge
  };
  struct Face {
    int SE; //starting edge
  };
public:
  //  mesh
  vector<Edge> _edges;
  vector<Vert> _verts;
  vector<Face> _faces;
  // points
  vector<Point3> _points;
public:
  TrMesh(const std::string& p): ComObject(p) {;}
  ~TrMesh() {;}
  //
  vector<Edge>& edges() { return _edges; }
  vector<Vert>& verts() { return _verts; }
  vector<Face>& faces() { return _faces; }
  vector<Point3>& points() { return _points; }
  //
  int setup(vector<Point3>& vv, vector<Index3>& fv);
  int numv() { return _verts.size(); }
  int nume() { return _edges.size(); }
  int numf() { return _faces.size(); }
  //  edge function
  bool Emajor(int E) {	 int CE = _edges[E].CE;	 return (E > CE);  }
  bool Einterior(int E) {	 int CE = _edges[E].CE;	 return (CE!=-1);  }
  int Econj(int E) {	 return _edges[E].CE;  }
  void E2Vs(int E, vector<int>& Vs) { //from E to get all Vs of E's face
    Vs.resize(3);//triangle
    int ME = E;          Vs[0] = _edges[ME].SV;
    ME = _edges[ME].NE;  Vs[1] = _edges[ME].SV;
    ME = _edges[ME].NE;  Vs[2] = _edges[ME].SV;
  }
  //  vert function
  bool Vclosed(int V) {	 int SE = _verts[V].SE;	 int CE = _edges[SE].CE;	 return (CE!=-1); } //for open ones, SE is the first one
  int Vvalence(int V) {
    int cnt = 0;	 int SE = _verts[V].SE;	 int ME = SE;	 bool start=true;
    while( (start==true || ME!=SE) && _edges[_edges[ME].PE].CE!=-1) {
      start = false;
      ME = _edges[_edges[ME].PE].CE;
      cnt ++;
    }
    return cnt;
  }
  //  face function
  void F2Vs(int F, vector<int>& Vs) {	 int E = _faces[F].SE;	 E2Vs(E, Vs);  }
  void F2Es(int F, vector<int>& Es) {
    Es.resize(3);
    int ME = _faces[F].SE;  Es[0] = ME;
    ME = _edges[ME].NE;     Es[1] = ME;
    ME = _edges[ME].NE;     Es[2] = ME;
  }
  //others
  void Fe2Fe(int F, int e, int& cF, int& ce) {
    vector<int> Es(3);	 F2Es(F, Es); //get all edges
    int E = Es[e]; //get eth edge
    int cE = _edges[E].CE; //find neighbor edge
    cF = _edges[cE].F; //find neighbor face
    F2Es(cF, Es); //get all edges
    for(int i=0; i<3; i++)		if(Es[i]==cE)		  ce = i; //compute ce
    return;
  }
  void Vf2Fv(int V, int f, int& F, int& v) {
    int ME = _verts[V].SE;
    for(int i=0; i<f; i++)
      ME = _edges[_edges[ME].PE].CE;
    F = _edges[ME].F; //face
    vector<int> Es(3);	 F2Es(F, Es);
    for(int i=0; i<3; i++)		if(Es[i]==ME)		  v = i;
    return;
  }
  void Fv2Vf(int F, int v, int& V, int& f) {
    vector<int> Vs(3);	 F2Vs(F, Vs);
    V = Vs[v];
    int ME = _verts[V].SE;
    int cnt = 0;
    while( _edges[ME].F != F) {
      ME = _edges[_edges[ME].PE].CE;
      cnt ++;
    }
    f = cnt;
    return;
  }
  //aux functions
  int nextvno(int vno) { return (vno+1)%3; }
  int prevvno(int vno) { return (vno+2)%3; }
  int nexteno(int eno) { return (eno+1)%3; }
  int preveno(int eno) { return (eno+2)%3; }
  int fromvno(int eno) { return eno; }
  int gotovno(int eno) { return (eno+1)%3; }
  //some useful functions
  int compute_dihedral( vector<double>& angles );
  int compute_interior( vector<double>& angles );
  int compute_area( vector<double>& areas );
};



#endif
