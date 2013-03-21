#include "wave3d.hpp"
#include "vecmatop.hpp"

using std::cerr;

//---------------------------------------------------------------------
int Wave3d::check(ParVec<int, cpx, PtPrtn>& den, ParVec<int, cpx, PtPrtn>& val, IntNumVec& chkkeys, double& relerr)
{
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  
  _self = this;
  int mpirank = this->mpirank();
  int mpisize = this->mpisize();
  
  ParVec<int, Point3, PtPrtn>& pos = (*_posptr);
  
  //1. get pos
  vector<int> all(1,1);
  vector<int> chkkeyvec;  for(int i=0; i<chkkeys.m(); i++)	chkkeyvec.push_back( chkkeys(i) );
  pos.getBegin(chkkeyvec, all);  pos.getEnd(all);
  //pos.get(chkkeyvec, all);
  vector<Point3> tmpsrcpos;
  for(map<int,Point3>::iterator mi=pos.lclmap().begin(); mi!=pos.lclmap().end(); mi++)
	if(pos.prtn().owner((*mi).first)==mpirank)
	  tmpsrcpos.push_back((*mi).second);
  vector<cpx> tmpsrcden;
  for(map<int,cpx>::iterator mi=den.lclmap().begin(); mi!=den.lclmap().end(); mi++)
	if(den.prtn().owner((*mi).first)==mpirank)
	  tmpsrcden.push_back((*mi).second);
  
  vector<Point3> tmptrgpos;
  for(int i=0; i<chkkeyvec.size(); i++)
	tmptrgpos.push_back( pos.access(chkkeyvec[i]) );
  
  DblNumMat srcpos(3, tmpsrcpos.size(), false, (double*)&(tmpsrcpos[0]));
  CpxNumVec srcden(tmpsrcden.size(), false, (cpx*)&(tmpsrcden[0]));
  DblNumMat trgpos(3, tmptrgpos.size(), false, (double*)&(tmptrgpos[0]));
  CpxNumVec trgval(tmptrgpos.size());
  
  CpxNumMat inter;
  iC( _knl.kernel(trgpos, srcpos, srcpos, inter) );
  iC( zgemv(1.0, inter, srcden, 0.0, trgval) );

  CpxNumVec allval(trgval.m());
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  iC( MPI_Allreduce(trgval.data(), allval.data(), trgval.m()*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) );
  
  //cerr<<trgval<<endl;
  //cerr<<allval<<endl;
  
  //2. get val
  val.getBegin(chkkeyvec, all);  val.getEnd(all);
  //val.get(chkkeyvec, all);
  CpxNumVec truval(chkkeyvec.size());
  for(int i=0; i<chkkeyvec.size(); i++)
	truval(i) = val.access(chkkeyvec[i]);
  //cerr<<truval<<endl;
  
  CpxNumVec errval(chkkeyvec.size());
  for(int i=0; i<chkkeyvec.size(); i++)
	errval(i) = allval(i) - truval(i);
  //cerr<<errval<<endl;
  
  double tn = sqrt( energy(truval) );
  double en = sqrt( energy(errval) );
  relerr = en/tn;
  
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  return 0;
}
