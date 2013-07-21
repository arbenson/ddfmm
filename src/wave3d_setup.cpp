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
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup");
#endif
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    _self = this;
    int mpirank = this->mpirank();
    //read optional data
    map<string,string>::iterator mi;
    mi = opts.find("-" + prefix() + "ACCU");
    if (mi != opts.end()) {
        istringstream ss(mi->second);
        ss >> _ACCU;
    }
    mi = opts.find("-" + prefix() + "NPQ");
    if (mi != opts.end()) {
        istringstream ss(mi->second);
        ss >> _NPQ;
    }
    mi = opts.find("-" + prefix() + "K");
    if (mi != opts.end()) {
        istringstream ss(mi->second);
        ss >> _K;
    }
    mi = opts.find("-" + prefix() + "ctr");
    if (mi != opts.end()) {
        istringstream ss(mi->second);
        double x, y, z;
        ss >> x >> y >> z;
        _ctr = Point3(x,y,z);
    }
    mi = opts.find("-" + prefix() + "ptsmax");
    if (mi != opts.end()) {
        istringstream ss(mi->second);
        ss >> _ptsmax;
    }
    mi = opts.find("-" + prefix() + "maxlevel");
    if (mi != opts.end()) {
        istringstream ss(mi->second);
        ss >> _maxlevel;
    }
    //
    if(mpirank == 0) {
        cerr << _K <<      " | "
             << _ACCU <<   " | "
             << _NPQ <<    " | "
             << _ctr <<    " | "
             << _ptsmax << " | "
             << _maxlevel
             << endl;
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
    _fplan = fftw_plan_dft_3d(2 * _P, 2 * _P, 2 * _P, (fftw_complex*) (_denfft.data()),
                              (fftw_complex*)(_denfft.data()), FFTW_FORWARD,
                              FFTW_MEASURE);
    iA(_fplan!=NULL);
    setvalue(_denfft,cpx(0,0));
    //
    _valfft.resize(2 * _P, 2 * _P, 2 * _P);
    _bplan = fftw_plan_dft_3d(2 * _P, 2 * _P, 2 * _P,
                              (fftw_complex*) (_valfft.data()),
                              (fftw_complex*)(_valfft.data()), FFTW_BACKWARD,
                              FFTW_ESTIMATE); 
    iA(_bplan!=NULL);
    setvalue(_valfft,cpx(0,0));
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_tree()
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_tree");
#endif
    int mpirank = this->mpirank();
    double eps = 1e-12;
    double K = this->K();
    // pos contains all data read by this processor
    ParVec<int, Point3, PtPrtn>& pos = (*_posptr);

    //1. call get
    // 1.  Get all of the geometry information needed for this processor
    vector<int> all(1,1);
    iC( pos.getBegin(&(Wave3d::setup_Q1_wrapper), all) );
    iC( pos.getEnd(all) );
    //iC( pos.get(&(Wave3d::setup_Q1_wrapper), all) );

    //2. _tree
    int numC = _geomprtn.m();
    int lvlC = cell_level();
    // generate cell level boxes, put them into queue
    Point3 bctr = ctr(); //OVERALL CENTER
    NumTns<BoxDat> cellboxtns(numC, numC, numC);
    for(map<int,Point3>::iterator mi=pos.lclmap().begin(); mi!=pos.lclmap().end(); mi++) {
        int key = mi->first;
        Point3 pos = mi->second;
        Index3 idx;
        for(int d = 0; d < 3; d++) {
            idx(d) = (int) floor(numC * ((pos(d) - bctr(d) + K / 2) / K));
            iA(idx(d) >= 0 && idx(d) < numC);
        }
        cellboxtns(idx(0),idx(1),idx(2)).ptidxvec().push_back( key ); //put the points in
    }
    //LEXING: NO MATTER WHETHER IT IS EMPTY OR NOT, ALWAYS IN
    queue< pair<BoxKey,BoxDat> > tmpq;
    for(int a = 0; a < numC; a++) {
        for(int b = 0; b < numC; b++) {
            for(int c = 0; c < numC; c++) {
                if(_geomprtn(a,b,c) == mpirank) {
                    BoxKey key(lvlC, Index3(a,b,c));
                    tmpq.push( pair<BoxKey,BoxDat>(key, cellboxtns(a,b,c)) );
                }
            }
        }
    }
    cellboxtns.resize(0,0,0);

    //-------tree, 
    while(!tmpq.empty()) {
        pair<BoxKey,BoxDat> curent = tmpq.front();
        tmpq.pop();
        BoxKey& curkey = curent.first;
        BoxDat& curdat = curent.second;
        //LEXING: VERY IMPORTANT
        if(curdat.ptidxvec().size() > 0) {
            curdat.tag() |=  WAVE3D_PTS;
        }
        bool action = (curkey.first <= unitlevel() && curdat.ptidxvec().size() > 0) ||
            (curdat.ptidxvec().size() > ptsmax() && curkey.first < maxlevel() - 1);
        if(action) {
            //1. subdivide to get new children
            NumTns<BoxDat> chdboxtns(2,2,2);
            Point3 curctr = center(curkey); //LEXING: VERY IMPORTANT
            for(int g = 0; g < curdat.ptidxvec().size(); g++) {
                int tmpidx = curdat.ptidxvec()[g];
                Point3 tmp = pos.access(tmpidx); //get position value
                Index3 idx;
                for(int d=0; d<3; d++) {
                    idx(d) = (tmp(d) >= curctr(d));
                }
                // put points to children
                chdboxtns(idx(0),idx(1),idx(2)).ptidxvec().push_back(tmpidx);
            }
            //2. put non-empty ones into queue
            for (int ind = 0; ind < NUM_DIRS; ind++) {
                int a = DIR_1(ind);
                int b = DIR_2(ind);
                int c = DIR_3(ind);
                //if(chdboxtns(a,b,c).ptidxvec().size()>0)
                BoxKey chdkey = this->chdkey(curkey, Index3(a,b,c));
                tmpq.push( pair<BoxKey,BoxDat>(chdkey, chdboxtns(a,b,c)) );
            }
            //4. clear my own ptidxvec vector
            curdat.ptidxvec().clear();
        } else {
            //1. copy data into _extpos
            curdat.extpos().resize(3, curdat.ptidxvec().size());
            for(int g = 0; g < curdat.ptidxvec().size(); g++) {
                int tmpidx = curdat.ptidxvec()[g];
                Point3 tmp = pos.access(tmpidx);
                for(int d=0; d<3; d++) {
                    curdat.extpos()(d,g) = tmp(d);
                }
            }
            //LEXING: VERY IMPORTANT
            curdat.tag() |= WAVE3D_TERMINAL;
        }
        //add my self into _tree
        _boxvec.insert(curkey, curdat); //LEXING: CHECK
    }

    //call get setup_Q2
    vector<int> mask1(BoxDat_Number,0);
    mask1[BoxDat_tag] = 1;
    iC( _boxvec.getBegin( &(Wave3d::setup_Q2_wrapper), mask1 ) );
    iC( _boxvec.getEnd( mask1 ) );
    //compute lists, low list and high list
    for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin();
        mi!=_boxvec.lclmap().end(); mi++) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if(own_box(curkey, mpirank) && has_pts(curdat)) { //LEXING: JUST COMPUTE MY OWN BOXES
            if(width(curkey) < 1-eps) { //LEXING: STRICTLY < 1
                // Low frequency regime
                iC( setup_tree_callowlist(curkey, curdat) );
            } else {
                // High frequency regime
                iC( setup_tree_calhghlist(curkey, curdat) );
            }
        }
    }

    //3. get extpos
    set<BoxKey> reqboxset;
    for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin();
        mi!=_boxvec.lclmap().end(); mi++) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if(has_pts(curdat) && own_box(curkey, mpirank)) {
            reqboxset.insert(curdat.undeidxvec().begin(), curdat.undeidxvec().end());
            reqboxset.insert(curdat.vndeidxvec().begin(), curdat.vndeidxvec().end());
            reqboxset.insert(curdat.wndeidxvec().begin(), curdat.wndeidxvec().end());
            reqboxset.insert(curdat.xndeidxvec().begin(), curdat.xndeidxvec().end());
        }
    }
    vector<BoxKey> reqbox;  reqbox.insert(reqbox.begin(), reqboxset.begin(), reqboxset.end());
    vector<int> mask2(BoxDat_Number,0);
    mask2[BoxDat_extpos] = 1;
    iC( _boxvec.getBegin(reqbox, mask2) );  iC( _boxvec.getEnd(mask2) );
    for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin();
        mi!=_boxvec.lclmap().end(); mi++) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (has_pts(curdat) && own_box(curkey, mpirank)) {
            for(vector<BoxKey>::iterator vi=curdat.vndeidxvec().begin();
                vi != curdat.vndeidxvec().end(); vi++) {
                BoxKey neikey = (*vi);
                BoxDat& neidat = _boxvec.access(neikey);
                neidat.fftnum() ++;
            }
        }
    }

    //4. dirupeqndenvec, dirdnchkvalvec
    //create
    for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin();
        mi!=_boxvec.lclmap().end(); mi++) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        double W = width(curkey);
        if(own_box(curkey, mpirank) && W>1-eps && has_pts(curdat)) {
            if(!iscell(curkey)) {
                BoxKey parkey = this->parkey(curkey);
                BoxDat& pardat = boxdata(parkey);
                for(set<Index3>::iterator si=pardat.outdirset().begin();
                    si != pardat.outdirset().end(); si++) {
                    Index3 nowdir = predir(*si);
                    curdat.outdirset().insert(nowdir);
                }
                for(set<Index3>::iterator si=pardat.incdirset().begin();
                    si != pardat.incdirset().end(); si++) {
                    Index3 nowdir = predir(*si);
                    curdat.incdirset().insert(nowdir);
                }
            }
            //go thrw
            Point3 curctr = center(curkey);
            for(map< Index3,vector<BoxKey> >::iterator mi=curdat.fndeidxvec().begin();
                mi!=curdat.fndeidxvec().end(); mi++) {
                vector<BoxKey>& tmplist = mi->second;
                for(int k = 0; k < tmplist.size(); k++) {
                    BoxKey othkey = tmplist[k];
                    Point3 othctr = center(othkey);
                    Point3 tmp = othctr - curctr;
                    tmp /= tmp.l2();
                    Index3 dir = nml2dir(tmp, W);
                    curdat.outdirset().insert(dir);

                    tmp = curctr - othctr;
                    tmp /= tmp.l2();
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
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_tree_callowlist");
#endif
    set<BoxKey> Uset, Vset, Wset, Xset;
    iA(!iscell(curkey)); //LEXING: DO NOT ALLOW WIDTH=1 BOXES TO BE CELLS
    Index3 curpth = curkey.second;
    BoxKey parkey = this->parkey(curkey); //the link to parent
    Index3 parpth = parkey.second;
    //
    int L = pow2(curkey.first);  //Index3 minpth(0,0,0);  //Index3 maxpth(L,L,L);
    int il, iu, jl,ju, kl, ku;
    il = max(2*parpth(0)-2,0);  iu = min(2*parpth(0)+4,L);
    jl = max(2*parpth(1)-2,0);  ju = min(2*parpth(1)+4,L);
    kl = max(2*parpth(2)-2,0);  ku = min(2*parpth(2)+4,L);
    for(int i = il; i < iu; i++) {
        for(int j = jl; j < ju; j++) {
            for(int k = kl; k < ku; k++) {
                Index3 trypth(i,j,k);
                if (trypth(0) == curpth(0)
		    && trypth(1) == curpth(1)
		    && trypth(2) == curpth(2)) {
                    continue;
                }
                BoxKey wntkey(curkey.first, trypth);
                //LEXING: LOOK FOR IT, DO NOT EXIST IF NO CELL BOX COVERING IT
                BoxKey reskey;
                bool found = setup_tree_find(wntkey, reskey);
                BoxDat& resdat = _boxvec.access(reskey);
                if (!found) {
                    continue;
                }
                bool adj = setup_tree_adjacent(reskey, curkey);

                if (reskey.first < curkey.first && has_pts(resdat)) {
                    if (!adj) {
                        Xset.insert(reskey);
                    } else if (isterminal(curdat)) {
                        Uset.insert(reskey);
                    }
                }

		if( reskey.first == curkey.first ) {
		    if(!adj) {
			Index3 bb = reskey.second - curkey.second;
			iA( bb.linfty()<=3 );
			if (has_pts(resdat)) {
			    Vset.insert(reskey);
			}
		    } else if(isterminal(curdat)) {
			queue<BoxKey> rest;
			rest.push(reskey);
			while(!rest.empty()) {
			    BoxKey fntkey = rest.front(); rest.pop();
			    BoxDat& fntdat = boxdata(fntkey);

			    bool adj = setup_tree_adjacent(fntkey, curkey);
			    if (!adj && has_pts(fntdat)) {
				Wset.insert(fntkey);
			    } 
			    if (adj && isterminal(fntdat) && has_pts(fntdat)) {
				Uset.insert(fntkey);
			    }
			    if (adj && !isterminal(fntdat)) {
				for (int ind = 0; ind < NUM_DIRS; ind++) {
				    rest.push( chdkey(fntkey, Index3(DIR_1(ind),
								     DIR_2(ind),
								     DIR_3(ind))) );
				}
			    }
			}
		    }
		}
	    }
	}
    }
    if (isterminal(curdat) && has_pts(curdat)) {
	Uset.insert(curkey);
    }
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
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_tree_calhghlist");
#endif
    Point3 curctr = center(curkey);
    double W = width(curkey);
    double eps = 1e-12;
    double D = W * W + W;
    double threshold = D - eps;
    if (iscell(curkey)) {
        //LEXING: CHECK THE FOLLOWING
        for(map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().begin(); iscell(mi->first); mi++) {
            BoxKey othkey = mi->first;
            BoxDat& othdat = _boxvec.access(othkey);
            if (has_pts(othdat)) {
                //LEXING: ALWAYS target - source //Point3 diff = curctr - center(othkey);
                Point3 diff = curctr - center(othkey);
                if(diff.l2() >= threshold) {
                    Index3 dir = nml2dir(diff / diff.l2(), W);
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
            for (int ind = 0; ind < NUM_DIRS; ind++) {
                BoxKey othkey = chdkey(trykey, Index3(DIR_1(ind), DIR_2(ind), DIR_3(ind)));
                BoxDat& othdat = _boxvec.access(othkey);
                if(!has_pts(othdat)) {
                    continue;
                }
                //LEXING: ALWAYS target - source
                Point3 diff = curctr - center(othkey);
                if(diff.l2() >= threshold) {
                    Index3 dir = nml2dir(diff/diff.l2(), W);
                    curdat.fndeidxvec()[dir].push_back(othkey);
                } else {
                    curdat.endeidxvec().push_back(othkey);
                }
            }
        }
    }
    return 0;
}

// ----------------------------------------------------------------------
bool Wave3d::setup_tree_find(BoxKey wntkey, BoxKey& trykey)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_tree_find");
#endif
    trykey = wntkey;
    while(!iscell(trykey)) {
        map<BoxKey,BoxDat>::iterator mi=_boxvec.lclmap().find(trykey);
        if(mi!=_boxvec.lclmap().end())
            return true; //found
        trykey = parkey(trykey);
    }
    map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().find(trykey);
    return (mi != _boxvec.lclmap().end());
}

// ----------------------------------------------------------------------
bool Wave3d::setup_tree_adjacent(BoxKey meekey, BoxKey youkey)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_tree_adjacent");
#endif
    int md = max(meekey.first,youkey.first);
    Index3 one(1,1,1);
    Index3 meectr(  (2 * meekey.second+one) * pow2(md - meekey.first)  );
    Index3 youctr(  (2 * youkey.second+one) * pow2(md - youkey.first)  );
    int meerad = pow2(md - meekey.first);
    int yourad = pow2(md - youkey.first);
    Index3 dif( ewabs(meectr-youctr) );
    int rad  = meerad + yourad;
    // return true iff at least one edge touch
    return dif[0] <= rad && dif[1] <= rad && dif[2] <= rad && dif.linfty() == rad;
}

//---------------------------------------------------------------------
int Wave3d::setup_Q1(int key, Point3& pos, vector<int>& pids)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_Q1");
#endif
    int numC = _geomprtn.m();
    Point3 center = ctr();
    Index3 idx;
    for(int d = 0; d < 3; d++) {
        idx(d) = (int) floor(numC * ((pos(d) - center(d) + _K / 2) / _K));
        iA(idx(d) >= 0 && idx(d) < numC);
    }
    pids.clear();
    pids.push_back( _geomprtn(idx(0),idx(1),idx(2)) ); //JUST ONE PROC
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_Q2(BoxKey boxkey, BoxDat& boxdat, vector<int>& pids)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_Q2");
#endif
    //for each ent, get all the pids that might need it
    int numC = _geomprtn.m();
    double widC = _K/numC;
    double W = width(boxkey);
    if (iscell(boxkey)) {
        //LEXING: CELL LEVEL BOXES ARE NEEDED FOR ALL CPUS
        pids.clear();
        for(int i = 0; i < mpisize(); i++) {
            pids.push_back(i);
        }
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
        for(int i = il; i < iu; i++) {
            for(int j = jl; j < ju; j++) {
                for(int k = kl; k < ku; k++) {
                    idset.insert( _geomprtn(i,j,k) );
                }
            }
        }
        pids.clear();
        pids.insert(pids.begin(), idset.begin(), idset.end());
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::setup_Q1_wrapper(int key, Point3& dat, vector<int>& pids)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_Q1_wrapper");
#endif
    return (Wave3d::_self)->setup_Q1(key, dat, pids);
}

int Wave3d::setup_Q2_wrapper(BoxKey key, BoxDat& dat, vector<int>& pids)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::setup_Q2_wrapper");
#endif
    return (Wave3d::_self)->setup_Q2(key, dat, pids);
}
