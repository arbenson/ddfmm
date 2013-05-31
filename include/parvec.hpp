#ifndef _PARVEC_HPP_
#define _PARVEC_HPP_

#include "commoninc.hpp"
#include "serialize.hpp"

using std::vector;
using std::map;
using std::pair;
using std::istream;
using std::ostream;
using std::istringstream;
using std::ostringstream;
using std::string;
using std::cerr;

//--------------------------------------------
template <class Key, class Data, class Partition>
class ParVec
{
public:
    map<Key,Data> _lclmap;
    Partition _prtn; //has function owner:Key->pid

    ParVec() {;}
    ~ParVec() {;}
    //
    map<Key,Data>& lclmap() { return _lclmap; }
    Partition& prtn() { return _prtn; }
    int insert(Key, Data&);
    
    // Get the data associated with Key from the local map
    Data& access(Key);

    // gather all entries st pid contains this proc
    // e2ps is a function that takes as input a Key-Data pair, and a vector of
    // ints to be filled with processor IDs.
    int getBegin(int (*e2ps)(Key, Data& ,vector<int>&), const vector<int>& mask);

    // gather all entries with key in keyvec
    int getBegin(vector<Key>& keyvec, const vector<int>& mask);
    
    int getEnd(const vector<int>& mask);
    
    //put data for all entries with key in keyvec
    int putBegin(vector<Key>& keyvec, const vector<int>& mask);    

    int putEnd(const vector<int>& mask);

    int expand( vector<Key>& keyvec); //allocate space for not-owned entries
    int discard(vector<Key>& keyvec); //remove non-owned entries

private:
    int mpirank() const { int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); return rank; }
    int mpisize() const { int size; MPI_Comm_size(MPI_COMM_WORLD, &size); return size; }
    int reset_vecs();
    int get_sizes(vector<int>& rszvec, vector<int>& sifvec);
    int make_buf_reqs(vector<int>& rszvec, vector<int>& sszvec);
    int strs2vec(vector<ostringstream*>& ossvec);

    //temporary data
    vector<int> _snbvec;
    vector<int> _rnbvec;
    vector< vector<char> > _sbufvec;
    vector< vector<char> > _rbufvec;
    MPI_Request *_reqs;
    MPI_Status  *_stats;
};

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::reset_vecs()
{
    int mpisize = this->mpisize();
    _snbvec.resize(mpisize, 0);
    for(int k = 0; k < mpisize; k++) {
        _snbvec[k] = 0;
    }
    _rnbvec.resize(mpisize, 0);
    for(int k = 0; k < mpisize; k++) {
        _rnbvec[k] = 0;
    }
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::get_sizes(vector<int>& rszvec, vector<int>& sszvec)
{
    int mpisize = this->mpisize();
    sszvec.resize(mpisize, 0);
    for(int k = 0; k < mpisize; k++) {
        sszvec[k] = _sbufvec[k].size();
    }
    vector<int> sifvec(2 * mpisize, 0);
    for(int k = 0; k < mpisize; k++) {
        sifvec[2 * k] = _snbvec[k];
        sifvec[2 * k + 1] = sszvec[k];
    }
    vector<int> rifvec(2 * mpisize, 0);
    iC( MPI_Alltoall( (void*)&(sifvec[0]), 2, MPI_INT, (void*)&(rifvec[0]), 2,
                      MPI_INT, MPI_COMM_WORLD ) );
    rszvec.resize(mpisize,0);
    for(int k = 0; k < mpisize; k++) {
        _rnbvec[k] = rifvec[2 * k];
        rszvec[k] = rifvec[2 * k + 1];
    }
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::make_buf_reqs(vector<int>& rszvec, vector<int>& sszvec)
{
    int mpisize = this->mpisize();
    for(int k = 0; k < mpisize; k++) {
        _rbufvec[k].resize(rszvec[k]);
    }
    for(int k = 0; k < mpisize; k++) {
        iC( MPI_Irecv( (void *)&(_rbufvec[k][0]), rszvec[k], MPI_BYTE, k, 0,
                       MPI_COMM_WORLD, &_reqs[2 * k] ) );
        iC( MPI_Isend( (void *)&(_sbufvec[k][0]), sszvec[k], MPI_BYTE, k, 0,
                       MPI_COMM_WORLD, &_reqs[2 * k + 1] ) );
    }
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::strs2vec(vector<ostringstream *>& ossvec)
{
    int mpisize = this->mpisize();
    for (int k = 0; k < mpisize; k++) {
        string tmp( ossvec[k]->str() );
        _sbufvec[k].clear();
        _sbufvec[k].insert(_sbufvec[k].end(), tmp.begin(), tmp.end());
    }
    for(int k = 0; k < mpisize; k++) {
        delete ossvec[k];
        ossvec[k] = NULL;
    }
    return 0;
}


//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::insert(Key key, Data& dat)
{
    int mpirank = this->mpirank();
    iA( _prtn.owner(key) == mpirank); //LEXING: VERY IMPORTANT
    _lclmap[key] = dat;
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
Data& ParVec<Key,Data,Partition>::access(Key key)
{
    typename map<Key,Data>::iterator mi=_lclmap.find(key);
    iA(mi != _lclmap.end());
    return mi->second;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::getBegin( int (*e2ps)(Key,Data&,vector<int>&),
                                          const vector<int>& mask )
{
    int mpirank = this->mpirank();
    int mpisize = this->mpisize();
    //---------
    reset_vecs();
    _sbufvec.resize(mpisize);
    _rbufvec.resize(mpisize);
    _reqs = (MPI_Request *) malloc(2*mpisize * sizeof(MPI_Request));
    _stats = (MPI_Status *) malloc(2*mpisize * sizeof(MPI_Status));
    //---------
    vector<ostringstream*> ossvec(mpisize);
    for(int k=0; k<mpisize; k++)        {
        ossvec[k] = new ostringstream();
    }

    //1. serialize
    for(typename map<Key,Data>::iterator mi=_lclmap.begin(); mi!=_lclmap.end(); mi++) {
        Key key = mi->first;
        const Data& dat = mi->second;
        if(_prtn.owner(key) == mpirank) {
            //ASK QUESTIONS
            vector<int> pids;
            int res = (*e2ps)(mi->first, mi->second, pids);
            for(int i = 0; i < pids.size(); i++) {
                int k = pids[i];
                if(k != mpirank) { //DO NOT SEND TO MYSELF
                    iC( serialize(key, *(ossvec[k]), mask) );
                    iC( serialize(dat, *(ossvec[k]), mask) );
                    _snbvec[k]++; //LEXING: VERY IMPORTANT
                }
            }
        }
    }

    // to vector
    strs2vec(ossvec);

    //2. all th sendsize of the message
    vector<int> sszvec;
    vector<int> rszvec;
    get_sizes(rszvec, sszvec);

    //3. allocate space, send and receive
    make_buf_reqs(rszvec, sszvec);
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::getBegin(vector<Key>& keyvec, const vector<int>& mask)
{
    int mpirank = this->mpirank();
    int mpisize = this->mpisize();
    //---------
    reset_vecs();
    _sbufvec.resize(mpisize);
    _rbufvec.resize(mpisize);
    _reqs = (MPI_Request *) malloc(2 * mpisize * sizeof(MPI_Request));
    _stats = (MPI_Status *) malloc(2 * mpisize * sizeof(MPI_Status));

    //1. go thrw the keyvec to partition them among other procs
    vector< vector<Key> > skeyvec(mpisize);
    for(int i = 0; i < keyvec.size(); i++) {
        Key key = keyvec[i];
        int owner = _prtn.owner(key);
        if(owner != mpirank) {
            skeyvec[owner].push_back(key);
        }
    }

    //2. setdn receive size of keyvec
    vector<int> sszvec(mpisize);
    vector<int> rszvec(mpisize);
    for(int k=0; k<mpisize; k++) {
        sszvec[k] = skeyvec[k].size();
    }
    iC( MPI_Alltoall( (void*)&(sszvec[0]), 1, MPI_INT, (void*)&(rszvec[0]), 1, MPI_INT, MPI_COMM_WORLD ) );

    //3. allocate space for the keys, send and receive
    vector< vector<Key> > rkeyvec(mpisize);
    for(int k=0; k<mpisize; k++) {
        rkeyvec[k].resize(rszvec[k]);
    }

    MPI_Request *reqs = (MPI_Request *) malloc(2*mpisize * sizeof(MPI_Request));
    MPI_Status  *stats = (MPI_Status *) malloc(2*mpisize * sizeof(MPI_Status));
    for(int k = 0; k < mpisize; k++) {
	iC( MPI_Irecv( (void*)&(rkeyvec[k][0]), rszvec[k]*sizeof(Key), MPI_BYTE,
                        k, 0, MPI_COMM_WORLD, &reqs[2 * k] ) );
	iC( MPI_Isend( (void*)&(skeyvec[k][0]), sszvec[k]*sizeof(Key), MPI_BYTE,
                       k, 0, MPI_COMM_WORLD, &reqs[2 * k + 1] ) );
    }
    iC( MPI_Waitall(2 * mpisize, &(reqs[0]), &(stats[0])) );
    free(reqs);
    free(stats);

    skeyvec.clear(); //save space

    //4. prepare the streams
    vector<ostringstream*> ossvec(mpisize);
    for(int k=0; k<mpisize; k++) {
        ossvec[k] = new ostringstream();
    }
    for(int k = 0; k < mpisize; k++) {
        for(int g = 0; g < rkeyvec[k].size(); g++) {
            Key curkey = rkeyvec[k][g];
            typename map<Key,Data>::iterator mi = _lclmap.find(curkey);
            iA( mi!=_lclmap.end() );
            iA( _prtn.owner(curkey) == mpirank );
            Key key = mi->first;
            const Data& dat = mi->second;
            iC( serialize(key, *(ossvec[k]), mask) );
            iC( serialize(dat, *(ossvec[k]), mask) );
            _snbvec[k]++; //LEXING: VERY IMPORTANT
        }
    }
    // to vector
    strs2vec(ossvec);

    //5. all the sendsize of the message
    get_sizes(rszvec, sszvec);

    //6. allocate space, send and receive
    make_buf_reqs(rszvec, sszvec);
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::getEnd( const vector<int>& mask )
{
    int mpisize = this->mpisize();
    //LEXING: SEPARATE HERE
    iC( MPI_Waitall(2*mpisize, &(_reqs[0]), &(_stats[0])) );
    free(_reqs);
    free(_stats);
    _sbufvec.clear(); //save space
    //4. write back
    //to stream
    vector<istringstream*> issvec(mpisize);
    for(int k=0; k<mpisize; k++) {
        issvec[k] = new istringstream();
    }
    for(int k=0; k<mpisize; k++) {
        string tmp(_rbufvec[k].begin(), _rbufvec[k].end());
        issvec[k]->str(tmp);
    }
    _rbufvec.clear(); //save space
  
    for (int k = 0; k < mpisize; k++) {
        for (int i = 0; i < _rnbvec[k]; i++) {
            Key key;  deserialize(key, *(issvec[k]), mask);
            typename map<Key,Data>::iterator mi=_lclmap.find(key);
            if (mi == _lclmap.end()) { //do not exist
                Data dat;
                deserialize(dat, *(issvec[k]), mask);
                _lclmap[key] = dat;
            } else { //exist already
                deserialize(mi->second, *(issvec[k]), mask);
            }
        }
    }
    for (int k = 0; k < mpisize; k++) {
        delete issvec[k];
        issvec[k] = NULL;
    }
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::putBegin(vector<Key>& keyvec, const vector<int>& mask)
{
    int mpirank = this->mpirank();
    int mpisize = this->mpisize();
    //---------
    reset_vecs();
    _sbufvec.resize(mpisize);
    _rbufvec.resize(mpisize);
    _reqs = (MPI_Request *) malloc(2*mpisize * sizeof(MPI_Request));
    _stats = (MPI_Status *) malloc(2*mpisize * sizeof(MPI_Status));
    //1.
    vector<ostringstream*> ossvec(mpisize);
    for (int k = 0; k < mpisize; k++) {
        ossvec[k] = new ostringstream();
    }

    //1. go thrw the keyvec to partition them among other procs
    for (int i=0; i<keyvec.size(); i++) {
        Key key = keyvec[i];
        int k = _prtn.owner(key); //the owner
        if (k != mpirank) {
            typename map<Key,Data>::iterator mi = _lclmap.find(key);
            iA( mi!=_lclmap.end() );
	    iA( key == mi->first );
            Data& dat = mi->second;
            iC( serialize(key, *(ossvec[k]), mask) );
            iC( serialize(dat, *(ossvec[k]), mask) );
            _snbvec[k]++;
        }
    }

    //2. to vector
    strs2vec(ossvec);

    //3. get size
    vector<int> sszvec;
    vector<int> rszvec;
    get_sizes(rszvec, sszvec);

    //4. allocate space, send and receive
    make_buf_reqs(rszvec, sszvec);

    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::putEnd( const vector<int>& mask )
{
    int mpirank = this->mpirank();
    int mpisize = this->mpisize();
    //LEXING: SEPARATE HERE
    iC( MPI_Waitall(2*mpisize, &(_reqs[0]), &(_stats[0])) );
    free(_reqs);
    free(_stats);
    _sbufvec.clear(); //save space
    //5. go thrw the messages and write back
    vector<istringstream*> issvec(mpisize);
    for(int k=0; k<mpisize; k++) {
        issvec[k] = new istringstream();
    }
    for(int k=0; k<mpisize; k++) {
        string tmp(_rbufvec[k].begin(), _rbufvec[k].end());
        issvec[k]->str(tmp);
    }
    _rbufvec.clear(); //save space
    for(int k=0; k<mpisize; k++) {
        for(int i=0; i<_rnbvec[k]; i++) {
            Key key;
            deserialize(key, *(issvec[k]), mask);
            iA( _prtn.owner(key) == mpirank );
            typename map<Key,Data>::iterator mi=_lclmap.find(key);
            iA( mi!=_lclmap.end() );
            deserialize(mi->second, *(issvec[k]), mask);
        }
    }
    for(int k=0; k<mpisize; k++) {
        delete issvec[k];
        issvec[k] = NULL;
    }
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::expand(vector<Key>& keyvec)
{
    Data dummy;
    for(int i=0; i<keyvec.size(); i++) {
        Key key = keyvec[i];
        typename map<Key,Data>::iterator mi=_lclmap.find(key);
        if(mi == _lclmap.end()) {
            _lclmap[key] = dummy;
        }
    }
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::discard(vector<Key>& keyvec)
{
    int mpirank = this->mpirank();
    for(int i=0; i<keyvec.size(); i++) {
        Key key = keyvec[i];
        if(_prtn.owner(key) != mpirank) {
            _lclmap.erase(key);
        }
    }
    return 0;
}

//-------------------
template<class Key, class Data, class Partition>
int serialize(const ParVec<Key,Data,Partition>& pv, ostream& os, const vector<int>& mask)
{
    serialize(pv._lclmap, os, mask);
    serialize(pv._prtn, os, mask);
    return 0;
}

template<class Key, class Data, class Partition>
int deserialize(ParVec<Key,Data,Partition>& pv, istream& is, const vector<int>& mask)
{
    deserialize(pv._lclmap, is, mask);
    deserialize(pv._prtn, is, mask);
    return 0;
}

#endif
