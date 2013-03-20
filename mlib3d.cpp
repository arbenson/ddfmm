#include "mlib3d.hpp"
#include "parallel.hpp"
#include "serialize.hpp"


using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::cerr;

//-----------------------------------
Mlib3d::Mlib3d(const string& p): ComObject(p)
{
}

//-----------------------------------
Mlib3d::~Mlib3d()
{
}

//-----------------------------------
int Mlib3d::setup(map<string,string>& opts)
{
  //get params
  map<string,string>::iterator mi;
  mi = opts.find("-" + prefix() + "NPQ");
  if(mi!=opts.end()) {
	istringstream ss((*mi).second);
	ss>>_NPQ;
  }
  mi = opts.find("-" + prefix() + "ldname");
  if(mi!=opts.end()) {
	istringstream ss((*mi).second);
	ss>>_ldname;
  }
  mi = opts.find("-" + prefix() + "hdname");
  if(mi!=opts.end()) {
	istringstream ss((*mi).second);
	ss>>_hdname;
  }
  
  //LEXING: read data in a shared way
  vector<int> all(1,1);
  istringstream liss;  iC( Shared_Read(_ldname, liss) );
  iC( deserialize(_w2ldmap, liss, all) );
  istringstream hiss;  iC( Shared_Read(_hdname, hiss) );
  iC( deserialize(_w2hdmap, hiss, all) );
  
  return 0;
}

//-----------------------------------
int Mlib3d::upward_lowfetch(double W, DblNumMat& uep, DblNumMat& ucp, NumVec<CpxNumMat>& uc2ue,
							NumTns<CpxNumMat>& ue2uc)
{
  { map<double,LowFreqEntry>::iterator mi = _w2ldmap.find(W);  iA(mi!=_w2ldmap.end());}
  LowFreqEntry& le = _w2ldmap[W];
  
  uep = le.uep();
  ucp = le.ucp();
  uc2ue.resize(3);
  uc2ue(0) = le.uc2ue()(0);
  uc2ue(1) = le.uc2ue()(1);
  uc2ue(2) = le.uc2ue()(2);
  
  DblNumMat uepchd;
  {
	LowFreqEntry& le = _w2ldmap[W/2];
	uepchd = le.uep();
  }
  ue2uc.resize(2,2,2);
  for(int a=0; a<2; a++)
	for(int b=0; b<2; b++)
	  for(int c=0; c<2; c++) {
		Point3 shift( (a-0.5)*W/2, (b-0.5)*W/2, (c-0.5)*W/2 );
		DblNumMat tmp(uepchd.m(), uepchd.n());
		for(int k=0; k<uepchd.n(); k++)		  for(int d=0; d<3; d++)			tmp(d,k) = uepchd(d,k) + shift(d);
		iC( _knl.kernel(ucp, tmp, tmp, ue2uc(a,b,c)) );
	  }
  
  return 0;
}

//-----------------------------------
int Mlib3d::dnward_lowfetch(double W, DblNumMat& dep, DblNumMat& dcp, NumVec<CpxNumMat>& dc2de,
							NumTns<CpxNumMat>& de2dc, NumTns<CpxNumTns>& ue2dc, DblNumMat& uep)
{
  { map<double,LowFreqEntry>::iterator mi = _w2ldmap.find(W);  iA(mi!=_w2ldmap.end());}
  LowFreqEntry& le = _w2ldmap[W];
  
  dep = le.ucp();
  dcp = le.uep();
  uep = le.uep();
  CpxNumMat tmp0( le.uc2ue()(0) );
  CpxNumMat tmp2( le.uc2ue()(2) );
  dc2de.resize(3);
  { CpxNumMat& trg = dc2de(0);  CpxNumMat& src = tmp2;  trg.resize(src.n(),src.m());  for(int i=0; i<trg.m(); i++)	for(int j=0; j<trg.n(); j++)	  trg(i,j) = src(j,i); }
  dc2de(1) = le.uc2ue()(1);
  { CpxNumMat& trg = dc2de(2);  CpxNumMat& src = tmp0;  trg.resize(src.n(),src.m());  for(int i=0; i<trg.m(); i++)	for(int j=0; j<trg.n(); j++)	  trg(i,j) = src(j,i); }
  DblNumMat dcpchd;
  {
	LowFreqEntry& le = _w2ldmap[W/2];
	dcpchd = le.uep();
  }
  
  de2dc.resize(2,2,2);
  for(int a=0; a<2; a++)
	for(int b=0; b<2; b++)
	  for(int c=0; c<2; c++) {
		Point3 shift( (a-0.5)*W/2, (b-0.5)*W/2, (c-0.5)*W/2 );
		DblNumMat tmp(dcpchd.m(), dcpchd.n());
		for(int k=0; k<dcpchd.n(); k++)		  for(int d=0; d<3; d++)			tmp(d,k) = dcpchd(d,k) + shift(d);
		iC( _knl.kernel(tmp, dep, dep, de2dc(a,b,c)) );
	  }
  
  ue2dc.resize(7,7,7);
  for(int a=0; a<7; a++)
	for(int b=0; b<7; b++)
	  for(int c=0; c<7; c++) {
		if(abs(a-3)>1 || abs(b-3)>1 || abs(c-3)>1)
		  ue2dc(a,b,c) = le.ue2dc()(a,b,c);
	  }
  return 0;
}

//-----------------------------------
int Mlib3d::upward_hghfetch(double W, Index3 dir, DblNumMat& uep, DblNumMat& ucp, NumVec<CpxNumMat>& uc2ue,
							NumTns<CpxNumMat>& ue2uc)
{
  { map< double, map<Index3,HghFreqDirEntry> >::iterator mi = _w2hdmap.find(W);  iA(mi!=_w2hdmap.end());}
  map<Index3,HghFreqDirEntry>& curmap = _w2hdmap[W];
  
  Index3 srt, sgn, prm;  iC( hghfetch_index3sort(dir, srt, sgn, prm) );
  { map<Index3,HghFreqDirEntry>::iterator mi = curmap.find(srt);  iA(mi!=curmap.end()); }
  HghFreqDirEntry& he = curmap[srt];
  DblNumMat ueptmp( he.uep() );
  DblNumMat ucptmp( he.ucp() );
  iC( hghfetch_shuffle(prm, sgn, ueptmp, uep) );
  iC( hghfetch_shuffle(prm, sgn, ucptmp, ucp) );
  uc2ue.resize(3);
  uc2ue(0) = he.uc2ue()(0); 
  uc2ue(1) = he.uc2ue()(1);
  uc2ue(2) = he.uc2ue()(2);
  
  DblNumMat uepchd;
  if(W==1.0) { //unit box	//map<double,int>::iterator mi = _w2ldmap.find(W/2);	iA(mi!=_w2ldmap.end());
	{ map<double,LowFreqEntry>::iterator mi = _w2ldmap.find(W/2);	iA(mi!=_w2ldmap.end()); }
	LowFreqEntry& le = _w2ldmap[W/2];
	uepchd = le.uep();
  } else { //large box
	{ map< double, map<Index3,HghFreqDirEntry> >::iterator mi = _w2hdmap.find(W/2);  iA(mi!=_w2hdmap.end());}
	map<Index3,HghFreqDirEntry>& curmap = _w2hdmap[W/2];
	
	Index3 pdr = predir(dir);	//cerr<<" 0 "<<dir<<" "<<pdr<<endl;
	Index3 srt, sgn, prm;  iC( hghfetch_index3sort(pdr, srt, sgn, prm) );
	{ map<Index3,HghFreqDirEntry>::iterator mi = curmap.find(srt);  iA(mi!=curmap.end()); }
	HghFreqDirEntry& he = curmap[srt];
	DblNumMat ueptmp( he.uep() );
	iC( hghfetch_shuffle(prm, sgn, ueptmp, uepchd) );
  }
  ue2uc.resize(2,2,2);
  for(int a=0; a<2; a++)
	for(int b=0; b<2; b++)
	  for(int c=0; c<2; c++) {
		Point3 shift( (a-0.5)*W/2, (b-0.5)*W/2, (c-0.5)*W/2 );
		DblNumMat tmp(uepchd.m(), uepchd.n());
		for(int k=0; k<uepchd.n(); k++)		  for(int d=0; d<3; d++)			tmp(d,k) = uepchd(d,k) + shift(d);
		iC( _knl.kernel(ucp, tmp, tmp, ue2uc(a,b,c)) );
	  }
  
  return 0;
}

//-----------------------------------
int Mlib3d::dnward_hghfetch(double W, Index3 dir, DblNumMat& dep, DblNumMat& dcp, NumVec<CpxNumMat>& dc2de,
							NumTns<CpxNumMat>& de2dc, DblNumMat& uep)
{
  { map< double, map<Index3,HghFreqDirEntry> >::iterator mi = _w2hdmap.find(W);  iA(mi!=_w2hdmap.end());}
  map<Index3,HghFreqDirEntry>& curmap = _w2hdmap[W];
  
  Index3 srt, sgn, prm;  iC( hghfetch_index3sort(dir, srt, sgn, prm) );
  { map<Index3,HghFreqDirEntry>::iterator mi = curmap.find(srt);  iA(mi!=curmap.end()); }
  HghFreqDirEntry& he = curmap[srt];
  
  DblNumMat ueptmp( he.uep() );
  DblNumMat ucptmp( he.ucp() );
  iC( hghfetch_shuffle(prm, sgn, ucptmp, dep) ); //ucp->dep
  for(int i=0; i<dep.n(); i++)	for(int d=0; d<3; d++)	  dep(d,i) = -dep(d,i); //negate
  iC( hghfetch_shuffle(prm, sgn, ueptmp, dcp) ); //uep->dcp
  for(int i=0; i<dcp.n(); i++)	for(int d=0; d<3; d++)	  dcp(d,i) = -dcp(d,i); //negate
  iC( hghfetch_shuffle(prm, sgn, ueptmp, uep) ); //uep->uep
  dc2de.resize(3);
  CpxNumMat tmp0( he.uc2ue()(0) );
  CpxNumMat tmp1( he.uc2ue()(1) );
  CpxNumMat tmp2( he.uc2ue()(2) );
  { CpxNumMat& trg = dc2de(0);  CpxNumMat& src = tmp2;  trg.resize(src.n(),src.m());  for(int i=0; i<trg.m(); i++)	for(int j=0; j<trg.n(); j++)	  trg(i,j) = src(j,i); }
  { CpxNumMat& trg = dc2de(1);  CpxNumMat& src = tmp1;  trg.resize(src.n(),src.m());  for(int i=0; i<trg.m(); i++)	for(int j=0; j<trg.n(); j++)	  trg(i,j) = src(j,i); }
  { CpxNumMat& trg = dc2de(2);  CpxNumMat& src = tmp0;  trg.resize(src.n(),src.m());  for(int i=0; i<trg.m(); i++)	for(int j=0; j<trg.n(); j++)	  trg(i,j) = src(j,i); }
  DblNumMat dcpchd;
  if(W==1.0) { //unit box
	{ map<double,LowFreqEntry>::iterator mi = _w2ldmap.find(W/2);	iA(mi!=_w2ldmap.end()); }
	LowFreqEntry& le = _w2ldmap[W/2];
	dcpchd = le.uep();
  } else { //large box
	{ map< double, map<Index3,HghFreqDirEntry> >::iterator mi = _w2hdmap.find(W/2);  iA(mi!=_w2hdmap.end());}
	map<Index3,HghFreqDirEntry>& curmap = _w2hdmap[W/2];
	
	Index3 pdr = predir(dir);
	Index3 srt, sgn, prm;  iC( hghfetch_index3sort(pdr, srt, sgn, prm) );
	{ map<Index3,HghFreqDirEntry>::iterator mi = curmap.find(srt);  iA(mi!=curmap.end()); }
	HghFreqDirEntry& he = curmap[srt];
	DblNumMat ueptmp( he.uep() );
	iC( hghfetch_shuffle(prm, sgn, ueptmp, dcpchd) );
	for(int i=0; i<dcpchd.n(); i++)	for(int d=0; d<3; d++)	  dcpchd(d,i) = -dcpchd(d,i); //negate
  }
  
  de2dc.resize(2,2,2);
  for(int a=0; a<2; a++)
	for(int b=0; b<2; b++)
	  for(int c=0; c<2; c++) {
		Point3 shift( (a-0.5)*W/2, (b-0.5)*W/2, (c-0.5)*W/2 );
		DblNumMat tmp(dcpchd.m(), dcpchd.n());
		for(int k=0; k<dcpchd.n(); k++)		  for(int d=0; d<3; d++)			tmp(d,k) = dcpchd(d,k) + shift(d);
		iC( _knl.kernel(tmp, dep, dep, de2dc(a,b,c)) );
	  }
  
  return 0;
}

//-----------------------------------
Index3 Mlib3d::predir(Index3 dir)
{
  int C = dir.linfty();
  int B = C/2;
  int midx = -1;
  for(int d=0; d<3; d++)
	if(abs(dir(d))==C) midx = d;
  //midx gives the direction
  Index3 res;
  for(int d=0; d<3; d++)	res(d) = (dir(d) + C - 1) / 2;
  for(int d=0; d<3; d++)	res(d) = 2*(res(d)/2) + 1 - B;
  res(midx) = dir(midx) / 2;
  return res;
}



//---------------------------------------------------------------------
int Mlib3d::hghfetch_shuffle(Index3 prm, Index3 sgn, DblNumMat& tmp, DblNumMat& res)
{
  res.resize(3, tmp.n());
  for(int k=0; k<res.n(); k++)
	for(int d=0; d<3; d++)
	  res(prm(d),k) = tmp(d,k);
  for(int k=0; k<res.n(); k++)
	for(int d=0; d<3; d++)
	  res(d,k) = sgn(d) * res(d,k);
  return 0;
}

//---------------------------------------------------------------------
int Mlib3d::hghfetch_index3sort(Index3 val, Index3& srt, Index3& sgn, Index3& prm)
{
  //make it positive
  for(int d=0; d<3; d++) {
	if(val(d)>=0)	  sgn(d) = 1;	else 	  sgn(d) = -1;
	val(d) = abs(val(d));
  }
  //sort
  int upp = val.linfty() + 1;
  for(int d=0; d<3; d++) {
	int minidx = d;
	int minval = val[d];
	for(int e=0; e<3; e++) {
	  if(val[e]<minval) {
		minval = val[e];
		minidx = e;
	  }
	}
	srt[d] = minval;
	prm[d] = minidx;
	val[minidx] = upp;
  }
  return 0;
}


//-------------------
int serialize(const LowFreqEntry& le, ostream& os, const vector<int>& mask)
{
  iC( serialize(le._uep, os, mask) );
  iC( serialize(le._ucp, os, mask) );
  iC( serialize(le._uc2ue, os, mask) );
  iC( serialize(le._ue2dc, os, mask) );
  return 0;
}
int deserialize(LowFreqEntry& le, istream& is, const vector<int>& mask)
{
  iC( deserialize(le._uep, is, mask) );
  iC( deserialize(le._ucp, is, mask) );
  iC( deserialize(le._uc2ue, is, mask) );
  iC( deserialize(le._ue2dc, is, mask) );
  return 0;
}

//-------------------
int serialize(const HghFreqDirEntry& he, ostream& os, const vector<int>& mask)
{
  iC( serialize(he._uep, os, mask) );
  iC( serialize(he._ucp, os, mask) );
  iC( serialize(he._uc2ue, os, mask) );
  return 0;
}

int deserialize(HghFreqDirEntry& he, istream& is, const vector<int>& mask)
{
  iC( deserialize(he._uep, is, mask) );
  iC( deserialize(he._ucp, is, mask) );
  iC( deserialize(he._uc2ue, is, mask) );
  return 0;
}
