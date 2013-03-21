#include "wave3d.hpp"
#include "parallel.hpp"

using namespace std;

int optionsCreate(int argc, char** argv, map<string,string>& options)
{
  //iA(argc%2==1);
  options.clear();
  //2. get extra data
  for(int k=1; k<argc; k=k+2) {
	options[ string(argv[k]) ] = string(argv[k+1]);
  }
  return 0;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //0. init and get options
  srand48(time(NULL));  iA(argc%1==0);
  clock_t ck0, ck1;
  time_t t0, t1;
  vector<int> all(1,1);
  map<string,string> opts;  optionsCreate(argc, argv, opts);
  map<string,string>::iterator mi;
  //1. read data
  ParVec<int, Point3, PtPrtn> pos;
  mi = opts.find("-posfile");  assert(mi!=opts.end());
  istringstream piss;  iC( Separate_Read((*mi).second, piss) );
  iC( deserialize(pos, piss, all) );
  vector<int>& tmpinfo = pos.prtn().ownerinfo();
  int numpts = tmpinfo[tmpinfo.size()-1]; //LEXING: THIS CONTAINS THE TOTAL NUMBER OF POINTS
  if(mpirank==0) {
    cout<<"Done reading pos "<<pos.lclmap().size()<<" "<<pos.prtn().ownerinfo().size()<<endl;
  }
  ParVec<int, cpx, PtPrtn> den;
  mi = opts.find("-denfile");  assert(mi!=opts.end());
  istringstream diss; iC( Separate_Read((*mi).second, diss) );
  iC( deserialize(den, diss, all) );
  if(mpirank==0) {
    cout<<"Done reading den "<<den.lclmap().size()<<" "<<den.prtn().ownerinfo().size()<<endl;
  }
  ParVec<int, cpx, PtPrtn> val; //preset val to be the same as den
  val = den;  //val.lclmap() = den.lclmap();  val.prtn() = den.prtn();
  if(mpirank==0) {
	cout<<"Done setting val "<<val.lclmap().size()<<" "<<val.prtn().ownerinfo().size()<<endl;
  }
  Kernel3d knl(KNL_HELM);
  mi = opts.find("-knl");
  if(mi!=opts.end()) {
	istringstream ss((*mi).second);
	int type;	ss>>type;
	knl.type() = type;
  }
  Mlib3d mlib("mlib3d_");
  mlib.knl() = knl;
  mlib.NPQ() = 4;
  iC( mlib.setup(opts) );
  if(mpirank==0) {
    cout<<"Done reading mlib "<<mlib.w2ldmap().size()<<" "<<mlib.w2hdmap().size()<<endl;
  }
  IntNumTns geomprtn;
  mi = opts.find("-geomprtn");  assert(mi!=opts.end());
  istringstream giss;  iC( Shared_Read((*mi).second, giss) );
  iC( deserialize(geomprtn, giss, all) );
  if(mpirank==0) {
    cout<<"Done reading geomprtn "<<geomprtn.m()<<" "<<geomprtn.n()<<" "<<geomprtn.p()<<endl;
    cout<<geomprtn<<endl;
  }
  Wave3d wave("wave3d_");
  wave.posptr() = &pos;
  wave.knl() = knl;
  wave.ACCU() = 1;
  wave.NPQ() = 4;
  wave.mlibptr() = &mlib;
  wave.geomprtn() = geomprtn;
  wave.K() = 64;
  wave.ctr() = Point3(0,0,0);
  wave.ptsmax() = 100; //LY: check
  //2. setup
  t0 = time(0);
  iC( wave.setup(opts) );
  t1 = time(0);
  if(mpirank==0) {	cout<<"wave setup used "<<difftime(t1,t0)<<"secs "<<endl;  }
  double time_setup = difftime(t1,t0);
  //3. eval
  double time_eval;
  for(int i=0; i<1; i++) {
	t0 = time(0);
	iC( wave.eval(den, val) );
	t1 = time(0);
	if(mpirank==0) {	cout<<"wave eval used "<<difftime(t1,t0)<<"secs "<<endl;  }
	time_eval = difftime(t1,t0);
  }
  ostringstream oss; iC( serialize(val, oss, all) );
  mi = opts.find("-valfile");  assert(mi!=opts.end());
  iC( Separate_Write((*mi).second, oss) );
  //4. check
  IntNumVec chkkeyvec;
  mi = opts.find("-chkfile");  assert(mi!=opts.end());
  istringstream iss;  iC( Shared_Read((*mi).second, iss) );
  iC( deserialize(chkkeyvec, iss, all) );
  int numchk = chkkeyvec.m();
  double relerr;
  ck0 = clock();
  iC( wave.check(den, val, chkkeyvec, relerr) );
  ck1 = clock();
  if(mpirank==0) {	 cout<<"wave check used "<<double(ck1-ck0)/CLOCKS_PER_SEC<<"secs "<<endl; }
  //
  double time_drct = (double(ck1-ck0)/CLOCKS_PER_SEC) * numpts / double(numchk);
  if(mpirank==0) {
	printf("----------------------\n");
	printf("RESULT\n");
	printf("K  %.2e\n", wave.K());
	printf("Ta %.2e\n", time_eval);
	printf("Td %.2e\n", time_drct);
	printf("Rt %.2e\n", time_drct/time_eval);
	printf("Ea %.2e\n", relerr);
	printf("----------------------\n");
  }
  //
  iC( MPI_Finalize() );
  return 0;
}
