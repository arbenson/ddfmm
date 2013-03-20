#include "wave3d.hpp"

using std::cerr;

//---------------------------------------------------------------------------
Wave3d* Wave3d::_self = NULL;

//-----------------------------------
Wave3d::Wave3d(const string& p): ComObject(p), _posptr(NULL), _mlibptr(NULL), _fplan(NULL), _bplan(NULL)
{
  _self = this;
}
//-----------------------------------
Wave3d::~Wave3d()
{
  _self = this;
  if(_fplan!=NULL) {	fftw_destroy_plan(_fplan);  }
  if(_bplan!=NULL) {	fftw_destroy_plan(_bplan);  }
}

//-----------------------------------
Index3 Wave3d::nml2dir(Point3 n, double W)
{
  int C = _NPQ * int(round(W));
  int B = C/2;
  int midx = 0;  double mval = abs(n(0));
  for(int d=1; d<3; d++)
	if(mval<abs(n(d))) {
	  midx = d;	  mval = abs(n(d));
	}
  //midx gives the direction (can be + or -)
  Point3 val = n / mval;
  for(int d=0; d<3; d++)	val(d) = atan(val(d));
  double Ang = (M_PI/2)/C;
  Index3 res;
  for(int d=0; d<3; d++) {
	int tmp = int(floor(val(d)/Ang));
	tmp = min(max(tmp,-B), B-1);
	res(d) = 2*tmp + 1;
  }
  res(midx) = C*int(round(val(midx))); //val(midx)==1 or -1
  return res;
}

//-----------------------------------
Index3 Wave3d::predir(Index3 dir)
{
  int C = dir.linfty();
  int B = C/2;
  int midx = -1;
  for(int d=0; d<3; d++)	if(abs(dir(d))==C) midx = d;
  assert(midx!=-1);
  //midx gives the direction
  Index3 res;
  for(int d=0; d<3; d++)	res(d) = (dir(d) + C - 1) / 2;
  for(int d=0; d<3; d++)	res(d) = 2*(res(d)/2) + 1 - B;
  res(midx) = dir(midx) / 2;
  return res;
}

//-----------------------------------
vector<Index3> Wave3d::chddir(Index3 dir)
{
  int C = dir.linfty();
  int midx = -1;
  vector<int> oidx;
  for(int d=0; d<3; d++)
	if(abs(dir(d))==C)
	  midx = d;
	else
	  oidx.push_back(d);
  vector<Index3> res;
  for(int a=0; a<2; a++)
	for(int b=0; b<2; b++) {
	  Index3 tmp = 2*dir;
	  tmp(oidx[0]) += 2*a-1;
	  tmp(oidx[1]) += 2*b-1;
	  res.push_back(tmp);
	}
  return res;
}

//-----------------------------------
double Wave3d::dir2width(Index3 dir)
{
  int C = dir.linfty();
  return double(C/_NPQ);
}

//--------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------
int serialize(const PtPrtn& val, ostream& os, const vector<int>& mask)
{
  serialize(val._ownerinfo, os, mask);
  return 0;
}

//-----------------------------------------------------------
int deserialize(PtPrtn& val, istream& is, const vector<int>& mask)
{
  deserialize(val._ownerinfo, is, mask);
  return 0;
}

//-----------------------------------------------------------
int serialize(const BoxDat& val, ostream& os, const vector<int>& mask)
{
  int i = 0;
  
  if(mask[i]==1) serialize(val._tag, os, mask);  i++;
  if(mask[i]==1) serialize(val._ptidxvec, os, mask);  i++;
  
  if(mask[i]==1) serialize(val._undeidxvec, os, mask);  i++;
  if(mask[i]==1) serialize(val._vndeidxvec, os, mask);  i++;
  if(mask[i]==1) serialize(val._wndeidxvec, os, mask);  i++;
  if(mask[i]==1) serialize(val._xndeidxvec, os, mask);  i++;
  if(mask[i]==1) serialize(val._endeidxvec, os, mask);  i++;
  if(mask[i]==1) serialize(val._fndeidxvec, os, mask);  i++;
  
  if(mask[i]==1) serialize(val._extpos, os, mask);  i++;
  if(mask[i]==1) serialize(val._extden, os, mask);  i++;
  if(mask[i]==1) serialize(val._upeqnden, os, mask);  i++;
  if(mask[i]==1) serialize(val._extval, os, mask);  i++;
  if(mask[i]==1) serialize(val._dnchkval, os, mask);  i++;
  
  if(mask[i]==1) serialize(val._upeqnden_fft, os, mask);  i++;
  if(mask[i]==1) serialize(val._incdirset, os, mask);  i++;
  if(mask[i]==1) serialize(val._outdirset, os, mask);  i++;
  if(mask[i]==1) serialize(val._fftnum, os, mask);  i++;
  if(mask[i]==1) serialize(val._fftcnt, os, mask);  i++;
  
  iA(i==BoxDat_Number);
  
  return 0;
}

//-----------------------------------------------------------
int deserialize(BoxDat& val, istream& is, const vector<int>& mask)
{
  int i = 0;
  
  if(mask[i]==1) deserialize(val._tag, is, mask);  i++;
  if(mask[i]==1) deserialize(val._ptidxvec, is, mask);  i++;

  if(mask[i]==1) deserialize(val._undeidxvec, is, mask);  i++;
  if(mask[i]==1) deserialize(val._vndeidxvec, is, mask);  i++;
  if(mask[i]==1) deserialize(val._wndeidxvec, is, mask);  i++;
  if(mask[i]==1) deserialize(val._xndeidxvec, is, mask);  i++;
  if(mask[i]==1) deserialize(val._endeidxvec, is, mask);  i++;
  if(mask[i]==1) deserialize(val._fndeidxvec, is, mask);  i++;

  if(mask[i]==1) deserialize(val._extpos, is, mask);  i++;
  if(mask[i]==1) deserialize(val._extden, is, mask);  i++;
  if(mask[i]==1) deserialize(val._upeqnden, is, mask);  i++;
  if(mask[i]==1) deserialize(val._extval, is, mask);  i++;
  if(mask[i]==1) deserialize(val._dnchkval, is, mask);  i++;
  
  if(mask[i]==1) deserialize(val._upeqnden_fft, is, mask);  i++;
  if(mask[i]==1) deserialize(val._incdirset, is, mask);  i++;
  if(mask[i]==1) deserialize(val._outdirset, is, mask);  i++;
  if(mask[i]==1) deserialize(val._fftnum, is, mask);  i++;
  if(mask[i]==1) deserialize(val._fftcnt, is, mask);  i++;
  
  iA(i==BoxDat_Number);
  
  return 0;
}

//-----------------------------------------------------------
int serialize(const BndDat& val, ostream& os, const vector<int>& mask)
{
  int i = 0;
  if(mask[i]==1) serialize(val._dirupeqnden, os, mask);  i++;
  if(mask[i]==1) serialize(val._dirdnchkval, os, mask);  i++;
  iA(i==BndDat_Number);
  return 0;
}
//-----------------------------------------------------------
int deserialize(BndDat& val, istream& is, const vector<int>& mask)
{
  int i = 0;
  if(mask[i]==1) deserialize(val._dirupeqnden, is, mask);  i++;
  if(mask[i]==1) deserialize(val._dirdnchkval, is, mask);  i++;
  iA(i==BndDat_Number);
  return 0;
}


