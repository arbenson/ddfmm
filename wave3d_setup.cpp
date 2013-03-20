#include "wave3d.hpp"

using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::set;
using std::queue;
using std::cerr;

//---------------------------------------------------------------------
int Wave3d::setup(map<string,string>& opts)
{
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  _self = this;
  int mpirank = this->mpirank();
  int mpisize = this->mpisize();
  //read optional data
  map<string,string>::iterator mi;
  mi = opts.find("-" + prefix() + "ACCU");
  if(mi!=opts.end()) {	istringstream ss((*mi).second);	ss>>_ACCU;  }
  mi = opts.find("-" + prefix() + "NPQ");
  if(mi!=opts.end()) {	istringstream ss((*mi).second);	ss>>_NPQ;  }
  mi = opts.find("-" + prefix() + "K");
  if(mi!=opts.end()) {	istringstream ss((*mi).second);	ss>>_K;  }
  mi = opts.find("-" + prefix() + "ctr");
  if(mi!=opts.end()) {	istringstream ss((*mi).second);	double x,y,z;	ss>>x>>y>>z;	_ctr = Point3(x,y,z);  }
  mi = opts.find("-" + prefix() + "ptsmax");
  if(mi!=opts.end()) {	istringstream ss((*mi).second);	ss>>_ptsmax;  }
  mi = opts.find("-" + prefix() + "maxlevel");
  if(mi!=opts.end()) {	istringstream ss((*mi).second);	ss>>_maxlevel;  }
  //
  if(mpirank==0) {
    cerr<<_K<<" | "<<_ACCU<<" | "<<_NPQ<<" | "<<_ctr<<" | "<<_ptsmax<<" | "<<_maxlevel<<endl;
  }
  //
  //create the parvecs
  BoxPrtn bp;  bp.ownerinfo() = _geomprtn;
  _boxvec.prtn() = bp;
  BndPrtn tp;  tp.ownerinfo() = _geomprtn;
  _bndvec.prtn() = tp;
  //generate octree
  iC( setup_tree() );
  //plans
  int _P = P();
  _denfft.resize(2*_P, 2*_P, 2*_P);
  _fplan = fftw_plan_dft_3d(2*_P,2*_P,2*_P,(fftw_complex*)(_denfft.data()),(fftw_complex*)(_denfft.data()), FFTW_FORWARD,FFTW_MEASURE);
  iA(_fplan!=NULL);
  setvalue(_denfft,cpx(0,0));
  //
  _valfft.resize(2*_P, 2*_P, 2*_P);
  _bplan = fftw_plan_dft_3d(2*_P,2*_P,2*_P,(fftw_complex*)(_valfft.data()),(fftw_complex*)(_valfft.data()), FFTW_BACKWARD,FFTW_ESTIMATE); 
  iA(_bplan!=NULL);
  setvalue(_valfft,cpx(0,0));
  return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_tree()
{
  int mpirank = this->mpirank();
  int mpisize = this->mpisize();
  double eps = 1e-12;
  double K = this->K();
  ParVec<int, Point3, PtPrtn>& pos = (*_posptr);
  //1. call get
  vector<int> all(1,1);
  iC( pos.getBegin(&(Wave3d::setup_Q1_wrapper), all) );  iC( pos.getEnd(all) );
  //iC( pos.get(&(Wave3d::setup_Q1_wrapper), all) );
  //2. _tree
  int numC = _geomprtn.m();
  int lvlC = celllevel();
  // generate cell level boxes, put them into queue
  Point3 bctr = ctr(); //OVERALL CENTER
  NumTns<BoxDat> cellboxtns(numC, numC, numC);
  for(map<int,Point3>::iterator mi=pos.lclmap().begin(); mi!=pos.lclmap().end(); mi++) {
    int key = (*mi).first;
    Point3 pos = (*mi).second;
    Index3 idx;
    for(int d=0; d<3; d++) {
      idx(d) = (int)floor(numC*((pos(d)-bctr(d)+K/2)/K));
      iA(idx(d)>=0 && idx(d)<numC);
    }
    cellboxtns(idx(0),idx(1),idx(2)).ptidxvec().push_back( key ); //put the points in
  }
  //LEXING: NO MATTER WHETHER IT IS EMPTY OR NOT, ALWAYS IN
  queue< pair<BoxKey,BoxDat> > tmpq;
  for(int a=0; a<numC; a++)
    for(int b=0; b<numC; b++)
      for(int c=0; c<numC; c++) {
	if(_geomprtn(a,b,c)==mpirank) {
	  BoxKey key;
	  key.first = lvlC;
	  key.second = Index3(a,b,c);
	  tmpq.push( pair<BoxKey,BoxDat>(key, cellboxtns(a,b,c)) );
	}
      }
  cellboxtns.resize(0,0,0);
  //-------tree, 
  while(tmpq.empty()==false) {
	pair<BoxKey,BoxDat> curent = tmpq.front();	tmpq.pop();
	BoxKey& curkey = curent.first;
	BoxDat& curdat = curent.second;
	//LEXING: VERY IMPORTANT
	if(curdat.ptidxvec().size()>0)
	  curdat.tag() = curdat.tag() | WAVE3D_PTS;
	bool action = (curkey.first<=unitlevel() && curdat.ptidxvec().size()>0) ||
	  (curdat.ptidxvec().size()>ptsmax() && curkey.first<maxlevel()-1);
	if(action) {
	  //1. subdivide to get new children
	  NumTns<BoxDat> chdboxtns(2,2,2);
	  Point3 curctr = center(curkey); //LEXING: VERY IMPORTANT
	  for(int g=0; g<curdat.ptidxvec().size(); g++) {
		int tmpidx = curdat.ptidxvec()[g];
		Point3 tmp = pos.access(tmpidx); //get position value
		Index3 idx;
		for(int d=0; d<3; d++)		  idx(d) = (tmp(d)>=curctr(d));
		chdboxtns(idx(0),idx(1),idx(2)).ptidxvec().push_back(tmpidx); //put points to children
	  }
	  //2. put non-empty ones into queue
	  for(int a=0; a<2; a++)
		for(int b=0; b<2; b++)
		  for(int c=0; c<2; c++) {
			//if(chdboxtns(a,b,c).ptidxvec().size()>0)
			BoxKey chdkey = this->chdkey(curkey, Index3(a,b,c));
			tmpq.push( pair<BoxKey,BoxDat>(chdkey, chdboxtns(a,b,c)) );
		  }
	  //4. clear my own ptidxvec vector
	  curdat.ptidxvec().clear();
	} else {
	  //1. copy data into _extpos
	  curdat.extpos().resize(3, curdat.ptidxvec().size());
	  for(int g=0; g<curdat.ptidxvec().size(); g++) {
	    int tmpidx = curdat.ptidxvec()[g];
	    Point3 tmp = pos.access(tmpidx);
	    for(int d=0; d<3; d++)		  curdat.extpos()(d,g) = tmp(d);
	  }
	  //LEXING: VERY IMPORTANT
	  curdat.tag() = curdat.tag() | WAVE3D_TERMINAL;
	}
	//add my self into _tree
	_boxvec.insert(curkey, curdat); //LEXING: CHECK
  }
  //call get setup_Q2
  vector<int> mask1(BoxDat_Number,0);
  mask1[BoxDat_tag] = 1;
  iC( _boxvec.getBegin( &(Wave3d::setup_Q2_wrapper), mask1 ) );  iC( _boxvec.getEnd( mask1 ) );
  //iC( _boxvec.get( &(Wave3d::setup_Q2_wrapper), mask1 ) );
  //compute lists, low list and high list
  for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin(); mi!=_boxvec.lclmap().end(); mi++) {
	BoxKey curkey = (*mi).first;
	BoxDat& curdat = (*mi).second;
	if(_boxvec.prtn().owner(curkey)==mpirank && (curdat.tag()&WAVE3D_PTS)) { //LEXING: JUST COMPUTE MY OWN BOXES
	  if(width(curkey)<1-eps) { //LEXING: STRICTLY < 1
		iC( setup_tree_callowlist(curkey, curdat) );
	  } else {
		iC( setup_tree_calhghlist(curkey, curdat) );
	  }
	}
  }
  //3. get extpos
  set<BoxKey> reqboxset;
  for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin(); mi!=_boxvec.lclmap().end(); mi++) {
    BoxKey curkey = (*mi).first;
    BoxDat& curdat = (*mi).second;
    if(curdat.tag() & WAVE3D_PTS) {
      if(_boxvec.prtn().owner(curkey)==mpirank) {
	for(int k=0; k<curdat.undeidxvec().size(); k++)		reqboxset.insert( curdat.undeidxvec()[k] );
	for(int k=0; k<curdat.vndeidxvec().size(); k++)		reqboxset.insert( curdat.vndeidxvec()[k] );
	for(int k=0; k<curdat.wndeidxvec().size(); k++)		reqboxset.insert( curdat.wndeidxvec()[k] );
	for(int k=0; k<curdat.xndeidxvec().size(); k++)		reqboxset.insert( curdat.xndeidxvec()[k] );
      }
    }
  }
  vector<BoxKey> reqbox;  reqbox.insert(reqbox.begin(), reqboxset.begin(), reqboxset.end());
  vector<int> mask2(BoxDat_Number,0);
  mask2[BoxDat_extpos] = 1;
  iC( _boxvec.getBegin(reqbox, mask2) );  iC( _boxvec.getEnd(mask2) );
  //iC( _boxvec.get(reqbox, mask2) );
  for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin(); mi!=_boxvec.lclmap().end(); mi++) {
	BoxKey curkey = (*mi).first;
	BoxDat& curdat = (*mi).second;
	if(curdat.tag() & WAVE3D_PTS) {
	  if(_boxvec.prtn().owner(curkey)==mpirank) {
		for(vector<BoxKey>::iterator vi=curdat.vndeidxvec().begin(); vi!=curdat.vndeidxvec().end(); vi++) {
		  BoxKey neikey = (*vi);
		  BoxDat& neidat = _boxvec.access(neikey);
		  neidat.fftnum() ++;
		}
	  }
	}
  }
  //4. dirupeqndenvec, dirdnchkvalvec
  //create
  for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin(); mi!=_boxvec.lclmap().end(); mi++) {
	BoxKey curkey = (*mi).first;
	BoxDat& curdat = (*mi).second;
	double W = width(curkey);
	if(_boxvec.prtn().owner(curkey)==mpirank && W>1-eps && (curdat.tag()&WAVE3D_PTS)) {
	  if(iscell(curkey)==false) {
		BoxKey parkey = this->parkey(curkey);
		BoxDat& pardat = boxdata(parkey);
		for(set<Index3>::iterator si=pardat.outdirset().begin(); si!=pardat.outdirset().end(); si++) {
		  Index3 nowdir = predir(*si);
		  curdat.outdirset().insert(nowdir);
		}
		for(set<Index3>::iterator si=pardat.incdirset().begin(); si!=pardat.incdirset().end(); si++) {
		  Index3 nowdir = predir(*si);
		  curdat.incdirset().insert(nowdir);
		}
	  }
	  //go thrw
	  Point3 curctr = center(curkey);
	  for(map< Index3,vector<BoxKey> >::iterator mi=curdat.fndeidxvec().begin(); mi!=curdat.fndeidxvec().end(); mi++) {
		vector<BoxKey>& tmplist = (*mi).second;
		for(int k=0; k<tmplist.size(); k++) {
		  BoxKey othkey = tmplist[k];
		  //BoxDat& othdat = _boxvec.access(othkey);
		  Point3 othctr = center(othkey);
		  Point3 tmp;
		  if(othkey>curkey)
			tmp = othctr-curctr;
		  else
			tmp = -(curctr-othctr);
		  tmp = tmp/tmp.l2();
		  Index3 dir = nml2dir(tmp, W);
		  curdat.outdirset().insert(dir);
		  if(curkey>othkey)
			tmp = curctr-othctr;
		  else
			tmp = -(othctr-curctr);
		  tmp = tmp/tmp.l2();
		  dir = nml2dir(tmp, W);
		  curdat.incdirset().insert(dir);
		}
	  }
	}
  }
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_tree_callowlist(BoxKey curkey, BoxDat& curdat)
{
  set<BoxKey> Uset, Vset, Wset, Xset;
  //const BoxKey& curkey = curent.first;
  //const BoxDat& curdat = curent.second;
  iA(iscell(curkey)==false); //LEXING: DO NOT ALLOW WIDTH=1 BOXES TO BE CELLS
  Index3 curpth = curkey.second;
  BoxKey parkey = this->parkey(curkey); //the link to parent
  Index3 parpth = parkey.second;
  //
  int L = pow2(curkey.first);  //Index3 minpth(0,0,0);  //Index3 maxpth(L,L,L);
  int il, iu, jl,ju, kl, ku;
  il = max(2*parpth(0)-2,0);	iu = min(2*parpth(0)+4,L);
  jl = max(2*parpth(1)-2,0);	ju = min(2*parpth(1)+4,L);
  kl = max(2*parpth(2)-2,0);	ku = min(2*parpth(2)+4,L);
  for(int i=il; i<iu; i++)
	for(int j=jl; j<ju; j++)
	  for(int k=kl; k<ku; k++) {
		Index3 trypth(i,j,k);
		if( (trypth(0)==curpth(0) && trypth(1)==curpth(1) && trypth(2)==curpth(2))==false ) {
		  BoxKey wntkey(curkey.first, trypth);
		  //LEXING: LOOK FOR IT, DO NOT EXIST IF NO CELL BOX COVERING IT
		  BoxKey reskey;
		  bool found = setup_tree_find(wntkey, reskey);
		  BoxDat& resdat = _boxvec.access(reskey);
		  if(found) {
			bool adj = setup_tree_adjacent(reskey, curkey);
			if( reskey.first<curkey.first ) {
			  if(adj) {
				if(isterminal(curdat)) {
				  if(resdat.tag()&WAVE3D_PTS)
					Uset.insert(reskey);
				}
			  } else {
				if(resdat.tag()&WAVE3D_PTS)
				  Xset.insert(reskey);
			  }
			}
			if( reskey.first==curkey.first ) {
			  if(!adj) {
				Index3 bb = reskey.second - curkey.second;				iA( bb.linfty()<=3 );
				if(resdat.tag()&WAVE3D_PTS)
				  Vset.insert(reskey);
			  } else {
				if(isterminal(curdat)) {
				  queue<BoxKey> rest;
				  rest.push(reskey);
				  while(rest.empty()==false) {
					BoxKey fntkey = rest.front(); rest.pop();
					BoxDat& fntdat = boxdata(fntkey);
					if(setup_tree_adjacent(fntkey, curkey)==false) {
					  if(fntdat.tag()&WAVE3D_PTS)
						Wset.insert(fntkey);
					} else {
					  if(isterminal(fntdat)) {
						if(fntdat.tag()&WAVE3D_PTS)
						  Uset.insert(fntkey);
					  } else {
						for(int a=0; a<2; a++)
						  for(int b=0; b<2; b++)
							for(int c=0; c<2; c++)
							  rest.push( chdkey(fntkey, Index3(a,b,c)) );
					  }
					}
				  }
				}
			  }
			}
		  }
		}
	  }
  if(isterminal(curdat))
	if(curdat.tag()&WAVE3D_PTS)
	  Uset.insert(curkey);
  //
  for(set<BoxKey>::iterator si=Uset.begin(); si!=Uset.end(); si++)
	curdat.undeidxvec().push_back(*si);
  for(set<BoxKey>::iterator si=Vset.begin(); si!=Vset.end(); si++)
	curdat.vndeidxvec().push_back(*si);
  for(set<BoxKey>::iterator si=Wset.begin(); si!=Wset.end(); si++)
	curdat.wndeidxvec().push_back(*si);
  for(set<BoxKey>::iterator si=Xset.begin(); si!=Xset.end(); si++)
	curdat.xndeidxvec().push_back(*si);
  return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_tree_calhghlist(BoxKey curkey, BoxDat& curdat)
{
  Point3 curctr = center(curkey);
  double W = width(curkey);
  int C = _NPQ*int(round(W));
  double eps = 1e-12;
  double D = W*W + W;
  if(iscell(curkey)==true) {
    //LEXING: CHECK THE FOLLOWING
    for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin(); iscell((*mi).first)==true; mi++) {
      BoxKey othkey = (*mi).first;
      BoxDat& othdat = _boxvec.access(othkey);
      if(othdat.tag() & WAVE3D_PTS) {
	//LEXING: ALWAYS target - source //Point3 diff = curctr - center(othkey);
	Point3 diff;
	if(curkey>othkey)
	  diff = curctr-center(othkey);
	else
	  diff = -(center(othkey)-curctr);
	if(diff.l2()>=D-eps) {
	  Index3 dir = nml2dir(diff/diff.l2(), W);
	  curdat.fndeidxvec()[dir].push_back(othkey);
	} else {
	  curdat.endeidxvec().push_back(othkey);
	}
      }
    }
  } else {
	BoxKey parkey = this->parkey(curkey);
	BoxDat& pardata = boxdata(parkey);
	for(int k=0; k<pardata.endeidxvec().size(); k++) {
	  BoxKey trykey = pardata.endeidxvec()[k];
	  BoxDat& trydata = boxdata(trykey);
	  for(int a=0; a<2; a++)
		for(int b=0; b<2; b++)
		  for(int c=0; c<2; c++) {
			BoxKey othkey = chdkey(trykey, Index3(a,b,c));
			BoxDat& othdat = _boxvec.access(othkey);
			if(othdat.tag() & WAVE3D_PTS) {
			  //LEXING: ALWAYS target - source
			  Point3 diff;
			  if(curkey>othkey)
				diff = curctr-center(othkey);
			  else
				diff = -(center(othkey)-curctr);
			  if(diff.l2()>=D-eps) {
				Index3 dir = nml2dir(diff/diff.l2(), W);
				curdat.fndeidxvec()[dir].push_back(othkey);
			  } else {
				curdat.endeidxvec().push_back(othkey);
			  }
			}
		  }
	}
  }
  return 0;
}

// ----------------------------------------------------------------------
bool Wave3d::setup_tree_find(BoxKey wntkey, BoxKey& trykey)
{
  trykey = wntkey;
  while(iscell(trykey)==false) {
	map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().find(trykey);
	if(mi!=_boxvec.lclmap().end())
	  return true; //found
	trykey = parkey(trykey);
  }
  map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().find(trykey);
  return (mi!=_boxvec.lclmap().end());
}

// ----------------------------------------------------------------------
bool Wave3d::setup_tree_adjacent(BoxKey meekey, BoxKey youkey)
{
  int md = max(meekey.first,youkey.first);
  Index3 one(1,1,1);
  Index3 meectr(  (2*meekey.second+one) * pow2(md - meekey.first)  );
  Index3 youctr(  (2*youkey.second+one) * pow2(md - youkey.first)  );
  int meerad = pow2(md - meekey.first);
  int yourad = pow2(md - youkey.first);
  Index3 dif( ewabs(meectr-youctr) );
  int rad  = meerad + yourad;
  return
	(dif[0]<=rad && dif[1]<=rad && dif[2]<=rad) &&
	( dif.linfty() == rad ); //at least one edge touch
}

//---------------------------------------------------------------------
int Wave3d::setup_Q1(int key, Point3& pos, vector<int>& pids)
{
  int numC = _geomprtn.m();
  Point3 center = ctr();
  Index3 idx;
  for(int d=0; d<3; d++) {
	idx(d) = (int)floor(numC*((pos(d)-center(d)+_K/2)/_K));
	iA(idx(d)>=0 && idx(d)<numC);
  }
  pids.clear();
  pids.push_back( _geomprtn(idx(0),idx(1),idx(2)) ); //JUST ONE PROC
  return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_Q2(BoxKey boxkey, BoxDat& boxdat, vector<int>& pids)
{
  //for each ent, get all the pids that might need it
  int numC = _geomprtn.m();
  double widC = _K/numC;
  double W = width(boxkey);
  if(iscell(boxkey)) {
	//LEXING: CELL LEVEL BOXES ARE NEEDED FOR ALL CPUS
	pids.clear();
	for(int i=0; i<mpisize(); i++)
	  pids.push_back(i);
  } else {
	set<int> idset;
	Point3 ctr = center(boxkey);
	double D = max(4*W*W+4*W, 1.0); //LEXING: THIS TAKE CARES THE LOW FREQUENCY PART
	int il = max((int)floor((ctr(0)+_K/2-D)/widC),0);
	int iu = min((int)ceil( (ctr(0)+_K/2+D)/widC),numC);
	int jl = max((int)floor((ctr(1)+_K/2-D)/widC),0);
	int ju = min((int)ceil( (ctr(1)+_K/2+D)/widC),numC);
	int kl = max((int)floor((ctr(2)+_K/2-D)/widC),0);
	int ku = min((int)ceil( (ctr(2)+_K/2+D)/widC),numC);
	//LEXING: IMPROVE THIS
	for(int i=il; i<iu; i++)
	  for(int j=jl; j<ju; j++)
		for(int k=kl; k<ku; k++)
		  idset.insert( _geomprtn(i,j,k) );
	pids.clear();
	pids.insert(pids.begin(), idset.begin(), idset.end());
  }
  return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_Q1_wrapper(int key, Point3& dat, vector<int>& pids)
{
  return (Wave3d::_self)->setup_Q1(key, dat, pids);
}

int Wave3d::setup_Q2_wrapper(BoxKey key, BoxDat& dat, vector<int>& pids)
{
  return (Wave3d::_self)->setup_Q2(key, dat, pids);
}

