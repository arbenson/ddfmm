#include "wave3d.hpp"
#include "vecmatop.hpp"

#include <algorithm>
#include <list>

using std::list;

#define DVMAX 400

ParData Wave3d::GatherParData(time_t t0, time_t t1) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::mean_var");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    double diff = difftime(t1, t0);
    double *rbuf = new double[mpisize];

    MPI_Gather((void *)&diff, 1, MPI_DOUBLE, rbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    ParData data = {0., 0., 0., 0.};
    if (mpirank == 0) {
	data.max = rbuf[0];
	data.min = rbuf[0];
        for (int i = 0; i < mpisize; i++) {
            data.mean += rbuf[i];
	    if (rbuf[i] > data.max)
		data.max = rbuf[i];
	    if (rbuf[i] < data.min)
		data.min = rbuf[i];
        }
        data.mean /= mpisize;

        for (int i = 0; i < mpisize; i++) {
            data.var += (rbuf[i] - data.mean) * (rbuf[i] - data.mean);
        }
        data.var /= (mpisize - 1);
    }

    delete[] rbuf;
    return data;
}

void PrintData(ParData data, std::string message) {
    int mpirank = getMPIRank();
    if (mpirank == 0) {
        cout << message << endl
	     << "mean: " << data.mean << endl
	     << "var: "  << data.var  << endl
	     << "max: "  << data.max  << endl
	     << "min: "  << data.min  << endl;
    }
}

bool CompareStackEntries(pair<double, Index3> a, pair<double, Index3> b) {
    return a.first < b.first;
}

int Wave3d::LowFreqUpwardPass(map< double, vector<BoxKey> >& ldmap,
                              set<BoxKey>& reqboxset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqUpwardPass");
#endif
    int mpirank = getMPIRank();
    if(mpirank == 0) {
        cout << "Beginning low frequency upward pass..." << endl;
    }

    time_t t0 = time(0);
    // For each box width in the low frequency regime that this processor
    // owns, evaluate upward.
    for(map< double, vector<BoxKey> >::iterator mi = ldmap.begin();
        mi != ldmap.end(); mi++) {
        iC( EvalUpwardLow(mi->first, mi->second, reqboxset) );
    }
    time_t t1 = time(0);
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    PrintData(GatherParData(t0, t1), "Low frequency upward pass");
    return 0;
}

int Wave3d::LowFreqDownwardComm(set<BoxKey>& reqboxset) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqDownwardComm");
#endif
    time_t t0 = time(0);
    vector<BoxKey> reqbox;
    reqbox.insert(reqbox.begin(), reqboxset.begin(), reqboxset.end());
    vector<int> mask(BoxDat_Number,0);
    mask[BoxDat_extden] = 1;
    mask[BoxDat_upeqnden] = 1;
    iC( _boxvec.getBegin(reqbox, mask) );
    iC( _boxvec.getEnd(mask) );
    time_t t1 = time(0);
    PrintData(GatherParData(t0, t1), "Low frequency downward communication");
    return 0;
}

int Wave3d::LowFreqDownwardPass(map< double, vector<BoxKey> >& ldmap) {
#ifndef RELEASE
    CallStackEntry entry("Wave3d::LowFreqDownwardPass");
#endif
    time_t t0 = time(0);
    for(map< double, vector<BoxKey> >::reverse_iterator mi = ldmap.rbegin();
        mi != ldmap.rend(); mi++) {
        iC( EvalDownwardLow(mi->first, mi->second) );
    }
    time_t t1 = time(0);
    PrintData(GatherParData(t0, t1), "Low frequency downward pass");

    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

#ifdef LIMITED_MEMORY
int Wave3d::HighFreqPass(map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap) {
# ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqPass");
# endif
    time_t t0, t1;
    int mpirank = getMPIRank();
    
    if(mpirank == 0) {
        cout << "Beginning high frequency pass..." << endl;
    }
    t0 = time(0);

    // Find all directions on the first level (width = 1)
    vector<Index3> basedirs;
    double local_max_W = 1;
    for (map<Index3, pair< vector<BoxKey>, vector<BoxKey> > >::iterator mi = hdmap.begin();
         mi != hdmap.end(); mi++) {
        Index3 dir = mi->first;
	double W = dir2width(dir);
        if (W == 1) {
            basedirs.push_back(dir);
        }
	if (W > local_max_W) {
	    local_max_W = W;
	}
    }

    double max_W = 1;
    MPI_Allreduce(&local_max_W, &max_W, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    set<BndKey> reqbndset;
    for(int i = 0; i < basedirs.size(); i++) {
	iC( EvalUpwardHighRecursive(1, basedirs[i], hdmap, reqbndset) );
    }
    t1 = time(0);
    PrintData(GatherParData(t0, t1), "High frequency upward pass");

    t0 = time(0);
    vector<BndKey> reqbnd;
    reqbnd.insert(reqbnd.begin(), reqbndset.begin(), reqbndset.end());
    vector<int> mask(BndDat_Number, 0);
    mask[BndDat_dirupeqnden] = 1;

    list< vector< pair<double, Index3> > > call_stacks;
    for (int i = 0; i < basedirs.size(); i++) {
	vector< pair<double, Index3> > curr_stack;
	Index3 dir = basedirs[i];
	iC( BuildDownwardHighCallStack(1, dir, hdmap, curr_stack) );
	std::sort(curr_stack.begin(), curr_stack.end(), &CompareStackEntries);
        call_stacks.push_back(curr_stack);
    }

    // Initial communication
    {
	vector<BndKey> curr_request;
	for (int i = 0; i < reqbnd.size(); i++) {
	    if (width(reqbnd[i].first) == max_W) {
		curr_request.push_back(reqbnd[i]);
	    }
	}
	_bndvec.initialize_data();
	iC( _bndvec.getBegin(curr_request, mask) );
	iC( _bndvec.getEnd(mask) );
        int recv = _bndvec.kbytes_received();
        int sent = _bndvec.kbytes_sent();
	int total_recv = 0;
	int total_sent = 0;
	MPI_Allreduce(&recv, &total_recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sent, &total_sent, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (mpirank == 0) {
	    cout << total_recv << " total kbytes received." << endl;
	    cout << total_sent << " total kbytes sent." << endl;
	}
    }

    for (double W = max_W; W >= 1; W /= 2) {
	if (W > 1) {
	    vector<BndKey> curr_request;
	    for (int i = 0; i < reqbnd.size(); i++) {
		if (width(reqbnd[i].first) == W / 2) {
		    curr_request.push_back(reqbnd[i]);
		}
	    }
	    _bndvec.initialize_data();
	    iC( _bndvec.getBegin(curr_request, mask) );
	    iC( _bndvec.getEnd(mask) );
	    int recv = _bndvec.kbytes_received();
	    int sent = _bndvec.kbytes_sent();
	    int total_recv = 0;
	    int total_sent = 0;
	    MPI_Allreduce(&recv, &total_recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(&sent, &total_sent, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	    if (mpirank == 0) {
		cout << total_recv << " total kbytes received." << endl;
		cout << total_sent << " total kbytes sent." << endl;
	    }
	}

	for (list< vector< pair<double, Index3> > >::iterator it = call_stacks.begin();
	     it != call_stacks.end(); ++it) {
	    while (!it->empty() && it->back().first == W) {
		Index3 dir = it->back().second;
                map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >::iterator mi = hdmap.find(dir);
		iA (mi != hdmap.end());
		EvalDownwardHigh(W, dir, mi->second);
                it->pop_back();
	    }
	}

    }
    t1 = time(0);

    PrintData(GatherParData(t0, t1), "High frequency downward pass");
    return 0;
}
#else
int Wave3d::HighFreqPass(map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap) {
# ifndef RELEASE
    CallStackEntry entry("Wave3d::HighFreqPass");
# endif
    time_t t0, t1;
    int mpirank = getMPIRank();
    
    if(mpirank == 0) {
        cout << "Beginning high frequency pass..." << endl;
    }
    t0 = time(0);

    // Find all directions on the first level (width = 1)
    vector<Index3> basedirs;
    for (map<Index3, pair< vector<BoxKey>, vector<BoxKey> > >::iterator mi = hdmap.begin();
         mi != hdmap.end(); mi++) {
        Index3 dir = mi->first;
        if (dir2width(dir) == 1) {
            basedirs.push_back(dir);
        }
    }

    set<BndKey> reqbndset;
    for(int i = 0; i < basedirs.size(); i++) {
	iC( EvalUpwardHighRecursive(1, basedirs[i], hdmap, reqbndset) );
    }
    t1 = time(0);
    PrintData(GatherParData(t0, t1), "High frequency upward pass");

    t0 = time(0);
    vector<int> mask(BndDat_Number,0);
    mask[BndDat_dirupeqnden] = 1;
    vector<BndKey> reqbnd;
    reqbnd.insert(reqbnd.begin(), reqbndset.begin(), reqbndset.end());
    iC( _bndvec.getBegin(reqbnd, mask) );
    iC( _bndvec.getEnd(mask) );
    t1 = time(0);
    PrintData(GatherParData(t0, t1), "High frequency communication");

    t0 = time(0);
    for (int i = 0; i < basedirs.size(); i++) {
	Index3 dir = basedirs[i]; // LEXING: PRE HERE
	iC( EvalDownwardHighRecursive(1, dir, hdmap) );
    }
    t1 = time(0);

    PrintData(GatherParData(t0, t1), "High frequency downward pass");
    return 0;
}
#endif

int Wave3d::GatherDensities(vector<int>& reqpts, ParVec<int,cpx,PtPrtn>& den) {
    int mpirank = getMPIRank();
    vector<int> all(1, 1);
    time_t t0 = time(0);
    iC( den.getBegin(reqpts, all) );
    iC( den.getEnd(all) );
    time_t t1 = time(0);
    if (mpirank == 0) {
        cout << "Density communication: " << difftime(t1, t0)
             << " secs" << endl;
    }
    return 0;
}

int Wave3d::ConstructMaps(map< double, vector<BoxKey> >& ldmap,
			  map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap) {
    int mpirank = getMPIRank();
    double eps = 1e-12;
    // construct maps, low frequency level by level, high frequency dir by dir
    // 
    // ldmap  maps box widths to a list of BoxKeys which correspond
    // to boxes in the low-frequency regime that are owned by this processor
    //
    // hdmap maps a direction to a pair of lists:
    //           1. List of BoxKeys of outgoing directions
    //           2. List of BoxKeys of incoming directions
    //
    // TODO (Austin): Make the data structure of hdmap clearer
    for (map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); mi++) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        double W = width(curkey);
        if (has_pts(curdat) && own_box(curkey, mpirank)) {
            // Boxes of width less than one that are nonempty and are owned
            // by this processor get put in the low-frequency map.
            if (W < 1 - eps) {
                ldmap[W].push_back(curkey);
            } else {
                // High frequency regime
                BndDat dummy;
                // For each outgoing direction of this box, add to the first list
                for (set<Index3>::iterator si = curdat.outdirset().begin();
                    si != curdat.outdirset().end(); si++) {
                    hdmap[*si].first.push_back(curkey);
                    // into bndvec
                    _bndvec.insert(BndKey(curkey, *si), dummy);
                }
		// For each incoming direction of this box, add to the second list
                for (set<Index3>::iterator si = curdat.incdirset().begin();
                    si != curdat.incdirset().end(); si++) {
                    hdmap[*si].second.push_back(curkey);
                    // into bndvec
                    _bndvec.insert(BndKey(curkey, *si), dummy);
                }
            }
        }
    }
    return 0;
}


//---------------------------------------------------------------------
int Wave3d::eval(ParVec<int,cpx,PtPrtn>& den, ParVec<int,cpx,PtPrtn>& val)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::eval");
#endif
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    _self = this;
    time_t t0, t1, t2, t3;
    int mpirank = getMPIRank();
    vector<int> all(1, 1);
    ParVec<int, Point3, PtPrtn>& pos = (*_posptr);
    // 1. Go through posptr to get nonlocal points
    vector<int> reqpts;
    for(map<int,Point3>::iterator mi = pos.lclmap().begin();
        mi != pos.lclmap().end(); mi++) {
        reqpts.push_back( mi->first );
    }

    GatherDensities(reqpts, den);

    // 3. compute extden using ptidxvec
    for(map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); mi++) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (has_pts(curdat) && own_box(curkey, mpirank) && isterminal(curdat)) {
            vector<int>& curpis = curdat.ptidxvec();
            CpxNumVec& extden = curdat.extden();
            extden.resize(curpis.size());
            for(int k = 0; k < curpis.size(); k++) {
                int poff = curpis[k];
                extden(k) = den.access(poff);
            }
        }
    }
    iC( den.discard(reqpts) );

    map< double, vector<BoxKey> > ldmap;
    map< Index3, pair< vector<BoxKey>, vector<BoxKey> > > hdmap;
    ConstructMaps(ldmap, hdmap);

    set<BoxKey> reqboxset;
    LowFreqUpwardPass(ldmap, reqboxset);

    HighFreqPass(hdmap);

    LowFreqDownwardComm(reqboxset);
    LowFreqDownwardPass(ldmap);

    //set val from extval
    vector<int> wrtpts;
    for(map<int,Point3>::iterator mi = pos.lclmap().begin();
        mi != pos.lclmap().end(); mi++) {
        if (pos.prtn().owner(mi->first) != mpirank) {
            wrtpts.push_back(mi->first);
        }
    }
    val.expand(wrtpts);
    for (map<BoxKey,BoxDat>::iterator mi = _boxvec.lclmap().begin();
        mi != _boxvec.lclmap().end(); mi++) {
        BoxKey curkey = mi->first;
        BoxDat& curdat = mi->second;
        if (has_pts(curdat) && own_box(curkey, mpirank) && isterminal(curdat)) {
            CpxNumVec& extval = curdat.extval();
            vector<int>& curpis = curdat.ptidxvec();
            for (int k = 0; k < curpis.size(); k++) {
                int poff = curpis[k];
                val.access(poff) = extval(k);
            }
        }
    }
    //call val->put
    val.putBegin(wrtpts, all);  val.putEnd(all);
    val.discard(wrtpts);
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalUpwardLow(double W, vector<BoxKey>& srcvec, set<BoxKey>& reqboxset)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalUpwardLow");
#endif
    DblNumMat uep;
    DblNumMat ucp;
    NumVec<CpxNumMat> uc2ue;
    NumTns<CpxNumMat> ue2uc;
    iC( _mlibptr->upwardLowFetch(W, uep, ucp, uc2ue, ue2uc) );
    //---------------
    int tdof = 1;
    for (int k = 0; k < srcvec.size(); k++) {
        BoxKey srckey = srcvec[k];
        BoxDat& srcdat = _boxvec.access(srckey);
        // If there are no points, continue to the next Box
        if (!has_pts(srcdat)) {
            continue;
        }
        //-----------------------------------------------------------------------------
        Point3 srcctr = center(srckey);
        //get array
        CpxNumVec upchkval(tdof*ucp.n());
        setvalue(upchkval,cpx(0,0));
        CpxNumVec& upeqnden = srcdat.upeqnden();
        //ue2dc
        if (isterminal(srcdat)) {
            DblNumMat upchkpos(ucp.m(), ucp.n());
            for (int k = 0; k < ucp.n(); k++) {
                for (int d = 0; d < dim(); d++) {
                    upchkpos(d,k) = ucp(d,k) + srcctr(d);
                }
            }
            //mul
            CpxNumMat mat;
            iC( _knl.kernel(upchkpos, srcdat.extpos(), srcdat.extpos(), mat) );
            iC( zgemv(1.0, mat, srcdat.extden(), 1.0, upchkval) );
        } else {
            for (int ind = 0; ind < NUM_DIRS; ind++) {
                int a = DIR_1(ind);
                int b = DIR_2(ind);
                int c = DIR_3(ind);
                BoxKey chdkey = this->chdkey(srckey, Index3(a, b, c));
                BoxDat& chddat = _boxvec.access(chdkey);
                if (has_pts(chddat)) {
                    iC( zgemv(1.0, ue2uc(a, b, c), chddat.upeqnden(), 1.0, upchkval) );
                }
            }
        }

        //uc2ue
        CpxNumMat& v  = uc2ue(0);
        CpxNumMat& is = uc2ue(1); //LEXING: it is stored as a matrix
        CpxNumMat& up = uc2ue(2);
        CpxNumVec mid(up.m());
        setvalue(mid,cpx(0,0));
        iC( zgemv(1.0, up, upchkval, 0.0, mid) );
        for (int k = 0; k < mid.m(); k++) {
            mid(k) = mid(k) * is(k,0);
        }
        upeqnden.resize(v.m());
        setvalue(upeqnden,cpx(0,0));
        iC( zgemv(1.0, v, mid, 0.0, upeqnden) );

        //-------------------------
        //EXTRA WORK, change role now
        // Add boxes in U, V, W, and X lists of trgdat to reqboxset.           
        BoxDat& trgdat = srcdat;
        reqboxset.insert(trgdat.undeidxvec().begin(), trgdat.undeidxvec().end());
        reqboxset.insert(trgdat.vndeidxvec().begin(), trgdat.vndeidxvec().end());
        reqboxset.insert(trgdat.wndeidxvec().begin(), trgdat.wndeidxvec().end());
        reqboxset.insert(trgdat.xndeidxvec().begin(), trgdat.xndeidxvec().end());
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalDownwardLow(double W, vector<BoxKey>& trgvec)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalDownwardLow");
#endif
    DblNumMat dep;
    DblNumMat dcp;
    NumVec<CpxNumMat> dc2de;
    NumTns<CpxNumMat> de2dc;
    NumTns<CpxNumTns> ue2dc;
    DblNumMat uep;
    iC( _mlibptr->downwardLowFetch(W, dep, dcp, dc2de, de2dc, ue2dc, uep) );
    //------------------
    int _P = P();
    for (int k = 0; k < trgvec.size(); k++) {
        BoxKey trgkey = trgvec[k];
        BoxDat& trgdat = _boxvec.access(trgkey);
        // If there are no points, continue to the next box.
        if (!has_pts(trgdat)) {
            continue;
        }
        //-----------------------------------------------------------------------------
        Point3 trgctr = center(trgkey);
        //array
        CpxNumVec& dnchkval = trgdat.dnchkval();
        if (dnchkval.m() == 0) {
            dnchkval.resize(dcp.n());
            setvalue(dnchkval,cpx(0,0));
        }
        if (trgdat.extval().m() == 0) {
            trgdat.extval().resize( trgdat.extpos().n() );
            setvalue(trgdat.extval(), cpx(0,0));
        }
        DblNumMat dnchkpos(dcp.m(), dcp.n());
        for (int k = 0; k < dcp.n(); k++) {
            for (int d = 0; d < dim(); d++) {
                dnchkpos(d,k) = dcp(d,k) + trgctr(d);
            }
        }
        // List computations
        iC( U_list_compute(trgdat) );
        iC( V_list_compute(trgdat, W, _P, trgctr, uep, dcp, dnchkval, ue2dc) );
        iC( W_list_compute(trgdat, W, uep) );
        iC( X_list_compute(trgdat, dcp, dnchkpos, dnchkval) );

        //-------------
        //dnchkval to dneqnden
        CpxNumMat& v  = dc2de(0);
        CpxNumMat& is = dc2de(1);
        CpxNumMat& up = dc2de(2);
        CpxNumVec mid(up.m());        setvalue(mid,cpx(0,0));
        iC( zgemv(1.0, up, dnchkval, 0.0, mid) );
        dnchkval.resize(0); //LEXING: SAVE SPACE
        for (int k = 0; k < mid.m(); k++) {
            mid(k) = mid(k) * is(k,0);
        }
        CpxNumVec dneqnden(v.m());
        setvalue(dneqnden,cpx(0,0));
        iC( zgemv(1.0, v, mid, 0.0, dneqnden) );
        //-------------
        //to children or to exact points
        if (isterminal(trgdat)) {
            DblNumMat dneqnpos(dep.m(), dep.n());
            for (int k = 0; k < dep.n(); k++) {
                for (int d = 0; d < dim(); d++) {
                    dneqnpos(d,k) = dep(d,k) + trgctr(d);
                }
            }
            //mul
            CpxNumMat mat;
            iC( _knl.kernel(trgdat.extpos(), dneqnpos, dneqnpos, mat) );
            iC( zgemv(1.0, mat, dneqnden, 1.0, trgdat.extval()) );
        } else {
            //put stuff to children
            for (int ind = 0; ind < NUM_DIRS; ind++) {
                int a = DIR_1(ind);
                int b = DIR_2(ind);
                int c = DIR_3(ind);
                BoxKey chdkey = this->chdkey(trgkey, Index3(a, b, c));
                BoxDat& chddat = _boxvec.access(chdkey);
                if (!has_pts(chddat)) {
                    continue;
                }
                //mul
                if (chddat.dnchkval().m() == 0) {
                    chddat.dnchkval().resize(de2dc(a,b,c).m());
                    setvalue(chddat.dnchkval(), cpx(0,0));
                }
                iC( zgemv(1.0, de2dc(a, b, c), dneqnden, 1.0, chddat.dnchkval()) );
            }
        }
    }
    return 0;
}

int Wave3d::U_list_compute(BoxDat& trgdat)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::U_list_compute");
#endif
    for(vector<BoxKey>::iterator vi = trgdat.undeidxvec().begin();
        vi != trgdat.undeidxvec().end(); vi++) {
        BoxKey neikey = (*vi);
        BoxDat& neidat = _boxvec.access(neikey);
        //mul
        CpxNumMat mat;
        iC( _knl.kernel(trgdat.extpos(), neidat.extpos(), neidat.extpos(), mat) );
        iC( zgemv(1.0, mat, neidat.extden(), 1.0, trgdat.extval()) );
    }
    return 0;
}

int Wave3d::V_list_compute(BoxDat& trgdat, double W, int _P, Point3& trgctr, DblNumMat& uep,
                           DblNumMat& dcp, CpxNumVec& dnchkval, NumTns<CpxNumTns>& ue2dc)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::V_list_compute");
#endif
    double step = W/(_P-1);
    setvalue(_valfft,cpx(0,0));
    //LEXING: SPECIAL
    for(vector<BoxKey>::iterator vi=trgdat.vndeidxvec().begin();
        vi!=trgdat.vndeidxvec().end(); vi++) {
        BoxKey neikey = (*vi);
        BoxDat& neidat = _boxvec.access(neikey);
        //mul
        Point3 neictr = center(neikey);         //double DD = neinde.width();
        Index3 idx;
        for(int d=0; d<dim(); d++) {
            idx(d) = int(round( (trgctr[d]-neictr[d])/W )); //LEXING:CHECK
        }
        //create if it is missing
        if(neidat.fftcnt()==0) {
            setvalue(_denfft, cpx(0,0));
            CpxNumVec& neiden = neidat.upeqnden();
            for(int k=0; k<uep.n(); k++) {
                int a = int( round((uep(0,k)+W/2)/step) ) + _P;
                int b = int( round((uep(1,k)+W/2)/step) ) + _P;
                int c = int( round((uep(2,k)+W/2)/step) ) + _P;
                _denfft(a,b,c) = neiden(k);
            }
            fftw_execute(_fplan);
            neidat.upeqnden_fft() = _denfft; //COPY to the right place
        }
        CpxNumTns& neidenfft = neidat.upeqnden_fft();
        //TODO: LEXING GET THE INTERACTION TENSOR
        CpxNumTns& inttns = ue2dc(idx[0]+3,idx[1]+3,idx[2]+3);
        for(int a = 0; a < 2 * _P; a++) {
            for(int b = 0; b < 2 * _P; b++) {
                for(int c = 0; c < 2 * _P; c++) {
                    _valfft(a,b,c) += (neidenfft(a,b,c)*inttns(a,b,c));
                }
            }
        }
        //clean if necessary
        neidat.fftcnt()++;
        if(neidat.fftcnt()==neidat.fftnum()) {
            neidat.upeqnden_fft().resize(0,0,0);
            neidat.fftcnt() = 0;//reset, LEXING
        }
    }
    fftw_execute(_bplan);
    //add back
    double coef = 1.0/(2*_P * 2*_P * 2*_P);
    for(int k=0; k<dcp.n(); k++) {
        int a = int( round((dcp(0,k)+W/2)/step) ) + _P;
        int b = int( round((dcp(1,k)+W/2)/step) ) + _P;
        int c = int( round((dcp(2,k)+W/2)/step) ) + _P;
        dnchkval(k) += (_valfft(a,b,c)*coef); //LEXING: VERY IMPORTANT
    }
    return 0;
}

int Wave3d::X_list_compute(BoxDat& trgdat, DblNumMat& dcp, DblNumMat& dnchkpos,
                           CpxNumVec& dnchkval)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::X_list_compute");
#endif
    for(vector<BoxKey>::iterator vi = trgdat.xndeidxvec().begin();
        vi != trgdat.xndeidxvec().end(); vi++) {
        BoxKey neikey = (*vi);
        BoxDat& neidat = _boxvec.access(neikey);
        Point3 neictr = center(neikey);
        if(isterminal(trgdat) && trgdat.extpos().n() < dcp.n()) {
            CpxNumMat mat;
            iC( _knl.kernel(trgdat.extpos(), neidat.extpos(), neidat.extpos(), mat) );
            iC( zgemv(1.0, mat, neidat.extden(), 1.0, trgdat.extval()) );
        } else {
            //mul
            CpxNumMat mat;
            iC( _knl.kernel(dnchkpos, neidat.extpos(), neidat.extpos(), mat) );
            iC( zgemv(1.0, mat, neidat.extden(), 1.0, dnchkval) );
        }
    }
    return 0;
}

int Wave3d::W_list_compute(BoxDat& trgdat, double W, DblNumMat& uep)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::W_list_compute");
#endif
    for(vector<BoxKey>::iterator vi = trgdat.wndeidxvec().begin();
        vi != trgdat.wndeidxvec().end(); vi++) {
        BoxKey neikey = (*vi);
        BoxDat& neidat = _boxvec.access(neikey);
        Point3 neictr = center(neikey);
        //upchkpos
        if (isterminal(neidat) && neidat.extpos().n()<uep.n()) {
            CpxNumMat mat;
            iC( _knl.kernel(trgdat.extpos(), neidat.extpos(), neidat.extpos(), mat) );
            iC( zgemv(1.0, mat, neidat.extden(), 1.0, trgdat.extval()) );
        } else {
            double coef = width(neikey) / W; //LEXING: SUPER IMPORTANT
            DblNumMat upeqnpos(uep.m(), uep.n()); //local version
            for(int k=0; k<uep.n(); k++) {
                for(int d=0; d<dim(); d++) {
                    upeqnpos(d,k) = coef*uep(d,k) + neictr(d);
                }
            }
            //mul
            CpxNumMat mat;
            iC( _knl.kernel(trgdat.extpos(), upeqnpos, upeqnpos, mat) );
            iC( zgemv(1.0, mat, neidat.upeqnden(), 1.0, trgdat.extval()) );
        }
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalUpwardHighRecursive(double W, Index3 nowdir,
        map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap,
        set<BndKey>& reqbndset)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalUpwardHighRecursive");
#endif
    map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >::iterator mi = hdmap.find(nowdir);
    if (mi != hdmap.end()) {
        iC( EvalUpwardHigh(W, nowdir, mi->second, reqbndset) );
        vector<Index3> dirvec = chddir(nowdir);
        for (int k = 0; k < dirvec.size(); k++) {
            iC( EvalUpwardHighRecursive(2 * W, dirvec[k], hdmap, reqbndset) );
        }
    }
    return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalDownwardHighRecursive(double W, Index3 nowdir,
        map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalDownwardHighRecursive");
#endif
    map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >::iterator mi = hdmap.find(nowdir);
    if (mi != hdmap.end()) {
        vector<Index3> dirvec = chddir(nowdir);
        for (int k = 0; k < dirvec.size(); k++) {
            iC( EvalDownwardHighRecursive(2 * W, dirvec[k], hdmap) );
        }
        iC( EvalDownwardHigh(W, nowdir, mi->second) );
    }
    return 0;
}

//---------------------------------------------------------------------
# ifdef LIMITED_MEMORY
int Wave3d::BuildDownwardHighCallStack(double W, Index3 nowdir,
				       map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >& hdmap,
				       vector< pair<double, Index3> >& call_stack)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::DownwardHighCallStack");
#endif
    map< Index3, pair< vector<BoxKey>, vector<BoxKey> > >::iterator mi = hdmap.find(nowdir);
    if (mi != hdmap.end()) {
        vector<Index3> dirvec = chddir(nowdir);
        for (int k = 0; k < dirvec.size(); k++) {
            iC( BuildDownwardHighCallStack(2 * W, dirvec[k], hdmap, call_stack) );
        }
	call_stack.push_back(pair<double, Index3>(W, nowdir));
    }
    return 0;
}
# endif

//---------------------------------------------------------------------
int Wave3d::EvalUpwardHigh(double W, Index3 dir,
        pair< vector<BoxKey>, vector<BoxKey> >& hdvecs, set<BndKey>& reqbndset)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalUpwardHigh");
#endif
    double eps = 1e-12;
    DblNumMat uep;
    DblNumMat ucp;
    NumVec<CpxNumMat> uc2ue;
    NumTns<CpxNumMat> ue2uc;
    iC( _mlibptr->upwardHighFetch(W, dir, uep, ucp, uc2ue, ue2uc) );
    //---------------
    vector<BoxKey>& srcvec = hdvecs.first;
    for(int k = 0; k < srcvec.size(); k++) {
        BoxKey srckey = srcvec[k];
        BoxDat& srcdat = _boxvec.access(srckey);
        // If there are no points, continue to the next box.
        if (!has_pts(srcdat)) {
            continue;
        }
        //-----------------------------------------------------------------------------
        Point3 srcctr = center(srckey);
        BndKey bndkey(srckey, dir);
        BndDat& bnddat = _bndvec.access( bndkey );
        CpxNumVec& upeqnden = bnddat.dirupeqnden();
        //eval
        CpxNumVec upchkval(ue2uc(0,0,0).m());
        setvalue(upchkval,cpx(0,0));
        // High-frequency M2M
        if(abs(W-1) < eps) {
            // The children boxes only have non-directional equivalent densities
            for (int ind = 0; ind < NUM_DIRS; ind++) {
                int a = DIR_1(ind);
                int b = DIR_2(ind);
                int c = DIR_3(ind);
                BoxKey chdkey = this->chdkey(srckey, Index3(a, b, c));
                BoxDat& chddat = _boxvec.access(chdkey);
                if(has_pts(chddat)) {
                    CpxNumVec& chdued = chddat.upeqnden();
                    iC( zgemv(1.0, ue2uc(a,b,c), chdued, 1.0, upchkval) );
                }
            }
        } else {
            Index3 pdir = predir(dir); // parent direction
            for (int ind = 0; ind < NUM_DIRS; ind++) {
                int a = DIR_1(ind);
                int b = DIR_2(ind);
                int c = DIR_3(ind);
                BoxKey chdkey = this->chdkey(srckey, Index3(a, b, c));
                BoxDat& chddat = _boxvec.access(chdkey);
                if(has_pts(chddat)) {
                    BndKey bndkey(chdkey, pdir);
                    BndDat& bnddat = _bndvec.access(bndkey);
                    CpxNumVec& chdued = bnddat.dirupeqnden();
                    iC( zgemv(1.0, ue2uc(a,b,c), chdued, 1.0, upchkval) );
                }
            }
        }

        // Upward check to upward equivalency (uc2ue)
        CpxNumMat& E1 = uc2ue(0);
        CpxNumMat& E2 = uc2ue(1);
        CpxNumMat& E3 = uc2ue(2);
        cpx dat0[DVMAX], dat1[DVMAX];
        CpxNumVec tmp0(E3.m(), false, dat0);
        iA(DVMAX>=E3.m());
        CpxNumVec tmp1(E2.m(), false, dat1);
        iA(DVMAX>=E2.m());
        upeqnden.resize(E1.m());
        setvalue(upeqnden,cpx(0,0));
        iC( zgemv(1.0, E3, upchkval, 0.0, tmp0) );
        iC( zgemv(1.0, E2, tmp0, 0.0, tmp1) );
        iC( zgemv(1.0, E1, tmp1, 0.0, upeqnden) );
    }

    iC( get_reqs(dir, hdvecs, reqbndset) );
    return 0;
}

int Wave3d::get_reqs(Index3 dir, pair< vector<BoxKey>, vector<BoxKey> >& hdvecs,
                     set<BndKey>& reqbndset)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::get_reqs");
#endif
  // Fill reqbndset
  vector<BoxKey>& trgvec = hdvecs.second;
  for(int k = 0; k < trgvec.size(); k++) {
      BoxKey trgkey = trgvec[k];
      BoxDat& trgdat = _boxvec.access(trgkey);
      vector<BoxKey>& tmpvec = trgdat.fndeidxvec()[dir];
      for(int i = 0; i < tmpvec.size(); i++) {
          BoxKey srckey = tmpvec[i];
          reqbndset.insert(BndKey(srckey, dir));
      }
  }
  return 0;
}

//---------------------------------------------------------------------
int Wave3d::EvalDownwardHigh(double W, Index3 dir,
			     pair< vector<BoxKey>, vector<BoxKey> >& hdvecs)
{
#ifndef RELEASE
    CallStackEntry entry("Wave3d::EvalDownwardHigh");
#endif
    int mpirank = getMPIRank();

    double eps = 1e-12;
    DblNumMat dep;
    DblNumMat dcp;
    NumVec<CpxNumMat> dc2de;
    NumTns<CpxNumMat> de2dc;
    DblNumMat uep;
    iC( _mlibptr->downwardHighFetch(W, dir, dep, dcp, dc2de, de2dc, uep) );
    //LEXING: IMPORTANT
    vector<BoxKey>& trgvec = hdvecs.second;

    for(int k = 0; k < trgvec.size(); k++) {
        BoxKey trgkey = trgvec[k];
        BoxDat& trgdat = _boxvec.access(trgkey);
        // If there are not points, continue to the next box.
        if (!has_pts(trgdat)) {
            continue;
        }
        //-----------------------------------------------------------------------------
        Point3 trgctr = center(trgkey);
        //1. mix
        //get target
        DblNumMat tmpdcp(dcp.m(),dcp.n());
        for (int k = 0; k < tmpdcp.n(); k++) {
            for (int d = 0; d < 3; d++) {
                tmpdcp(d,k) = dcp(d,k) + trgctr(d);
            }
        }
        BndKey bndkey(trgkey, dir);
        BndDat& bnddat = _bndvec.access(bndkey);
        CpxNumVec& dcv = bnddat.dirdnchkval();
        vector<BoxKey>& tmpvec = trgdat.fndeidxvec()[dir];
        for (int i = 0; i < tmpvec.size(); i++) {
            BoxKey srckey = tmpvec[i];
            Point3 srcctr = center(srckey);
            //difference vector
            Point3 diff = trgctr - srcctr;
            diff /= diff.l2(); //LEXING: see wave3d_setup.cpp
            iA( nml2dir(diff, W) == dir );
            //get source
            DblNumMat tmpuep(uep.m(),uep.n());
            for (int k = 0; k < tmpuep.n(); k++) {
                for (int d = 0; d < 3; d++) {
                    tmpuep(d,k) = uep(d,k) + srcctr(d);
                }
            }
            BndKey bndkey(srckey, dir);
            BndDat& bnddat = _bndvec.access(bndkey);
            CpxNumVec& ued = bnddat.dirupeqnden();
            //mateix
            CpxNumMat Mts;
	    iC( _knl.kernel(tmpdcp, tmpuep, tmpuep, Mts) );
            //allocate space if necessary
            if (dcv.m() == 0) {
                dcv.resize(Mts.m());
                setvalue(dcv,cpx(0,0)); //LEXING: CHECK
            }
	    //iC( ued.m() != 0 );
	    if (ued.m() == 0) {
		ued.resize(Mts.n());
		setvalue(ued, cpx(0, 0));
	    }
            iC( zgemv(1.0, Mts, ued, 1.0, dcv) );
        }
        //2. to children
        CpxNumVec& dnchkval = dcv;
        //dc2de
        CpxNumMat& E1 = dc2de(0);
        CpxNumMat& E2 = dc2de(1);
        CpxNumMat& E3 = dc2de(2);
        cpx dat0[DVMAX], dat1[DVMAX], dat2[DVMAX];
        CpxNumVec tmp0(E3.m(), false, dat0);
        CpxNumVec tmp1(E2.m(), false, dat1);
        CpxNumVec dneqnden(E1.m(), false, dat2);
        iC( zgemv(1.0, E3, dnchkval, 0.0, tmp0) );
        iC( zgemv(1.0, E2, tmp0, 0.0, tmp1) );
        iC( zgemv(1.0, E1, tmp1, 0.0, dneqnden) );
        dnchkval.resize(0); //LEXING: SAVE SPACE
        //eval
        if (abs(W - 1)<eps) {
            for (int ind = 0; ind < NUM_DIRS; ind++) {
                int a = DIR_1(ind);
                int b = DIR_2(ind);
                int c = DIR_3(ind);             
                BoxKey chdkey = this->chdkey(trgkey, Index3(a,b,c));
                BoxDat& chddat = _boxvec.access(chdkey);
                if (!has_pts(chddat)) {
                    continue;
                }
                CpxNumVec& chddcv = chddat.dnchkval();
                if (chddcv.m() == 0) {
                    chddcv.resize(de2dc(a,b,c).m());
                    setvalue(chddcv,cpx(0,0));
                }
                iC( zgemv(1.0, de2dc(a,b,c), dneqnden, 1.0, chddcv) );
            }
        } else {
            Index3 pdir = predir(dir); //LEXING: CHECK
            for (int ind = 0; ind < NUM_DIRS; ind++) {
                int a = DIR_1(ind);
                int b = DIR_2(ind);
                int c = DIR_3(ind);             
                BoxKey chdkey = this->chdkey(trgkey, Index3(a,b,c));
                BoxDat& chddat = _boxvec.access(chdkey);
                if (!has_pts(chddat)) {
                    continue;
                }
                BndKey bndkey(chdkey, pdir);
                BndDat& bnddat = _bndvec.access(bndkey);
                CpxNumVec& chddcv = bnddat.dirdnchkval();
                if(chddcv.m()==0) {
                    chddcv.resize(de2dc(a,b,c).m());
                    setvalue(chddcv,cpx(0,0));
                }
                iC( zgemv(1.0, de2dc(a,b,c), dneqnden, 1.0, chddcv) );
            }
        }
    }
    //-----------------
    //EXTRA WORK, change role
    vector<BoxKey>& srcvec = hdvecs.first;
    for (int k = 0; k < srcvec.size(); k++) {
        BoxKey srckey = srcvec[k];
        BoxDat& srcdat = _boxvec.access(srckey);
        if (!has_pts(srcdat)) {
            continue;
        }
        BndKey bndkey(srckey, dir);
        BndDat& bnddat = _bndvec.access( bndkey );
        bnddat.dirupeqnden().resize(0);
    }
    for (int k = 0; k < trgvec.size(); k++) {
        BoxKey trgkey = trgvec[k];
        BoxDat& trgdat = _boxvec.access(trgkey);
        if (!has_pts(trgdat)) {
            continue;
        }
        vector<BoxKey>& tmpvec = trgdat.fndeidxvec()[dir];
        for (int i = 0; i < tmpvec.size(); i++) {
            BoxKey srckey = tmpvec[i];
            BndKey bndkey(srckey, dir);
            BndDat& bnddat = _bndvec.access(bndkey);
            bnddat.dirupeqnden().resize(0);
        }
    }
    return 0;
}
