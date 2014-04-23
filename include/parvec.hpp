/* Distributed Directional Fast Multipole Method
   Copyright (C) 2014 Austin Benson, Lexing Ying, and Jack Poulson

 This file is part of DDFMM.

    DDFMM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DDFMM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DDFMM.  If not, see <http://www.gnu.org/licenses/>. */
#ifndef _PARVEC_HPP_
#define _PARVEC_HPP_

#include "commoninc.hpp"
#include "serialize.hpp"

//--------------------------------------------
template <class Key, class Data, class Partition>
class ParVec {
public:
    std::map<Key,Data> _lclmap;
    Partition _prtn; // has function owner:Key->pid

    ParVec() {;}
    ~ParVec() {;}
    //
    std::map<Key,Data>& lclmap() { return _lclmap; }
    Partition& prtn() { return _prtn; }
    int insert(Key, Data&);
    
    // Get the data associated with Key from the local map
    Data& access(Key);

    std::pair<bool, Data&> contains(Key);

    // gather all entries st pid contains this proc
    // e2ps is a function that takes as input a Key-Data pair, and a vector of
    // ints to be filled with processor IDs.
    int getBegin(int (*e2ps)(Key, Data&, std::vector<int>&), const std::vector<int>& mask);

    // gather all entries with key in keyvec
    int getBegin(std::vector<Key>& keyvec, const std::vector<int>& mask);
    
    int getEnd(const std::vector<int>& mask);
    
    // put data for all entries with key in keyvec
    int putBegin(std::vector<Key>& keyvec, const std::vector<int>& mask);    

    int putEnd(const std::vector<int>& mask);

    // Allocate space for non-owned data for keys in keyvec.
    int expand(std::vector<Key>& keyvec);
    // Remove non-owned data for keys in keyvec.
    int discard(std::vector<Key>& keyvec);

    // Functions for keeping track of the amount of data that is sent
    // and received.
    int initialize_data() {
        _kbytes_received = 0;
        _kbytes_sent = 0;
    }
    int kbytes_received() { return _kbytes_received; }
    int kbytes_sent() { return _kbytes_sent; }

private:
    int resetVecs();
    int getSizes(std::vector<int>& rszvec, std::vector<int>& sifvec);
    int makeBufReqs(std::vector<int>& rszvec, std::vector<int>& sszvec);
    int strs2vec(std::vector<std::ostringstream*>& ossvec);

    // temporary data
    std::vector<int> _snbvec;
    std::vector<int> _rnbvec;
    std::vector< std::vector<char> > _sbufvec;
    std::vector< std::vector<char> > _rbufvec;
    MPI_Request *_reqs;
    MPI_Status  *_stats;
    std::string _tag;

    int _kbytes_received;
    int _kbytes_sent;
};

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::resetVecs() {
#ifndef RELEASE
    CallStackEntry entry("ParVec::resetVecs");
#endif
    int mpisize = getMPISize();
    _snbvec.resize(mpisize, 0);
    for (int k = 0; k < mpisize; ++k) {
        _snbvec[k] = 0;
    }
    _rnbvec.resize(mpisize, 0);
    for (int k = 0; k < mpisize; ++k) {
        _rnbvec[k] = 0;
    }
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::getSizes(std::vector<int>& rszvec,
                                         std::vector<int>& sszvec) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::getSizes");
#endif
    int mpisize = getMPISize();
    sszvec.resize(mpisize, 0);
    for (int k = 0; k < mpisize; ++k) {
        sszvec[k] = _sbufvec[k].size();
    }
    std::vector<int> sifvec(2 * mpisize, 0);
    for (int k = 0; k < mpisize; ++k) {
        sifvec[2 * k] = _snbvec[k];
        sifvec[2 * k + 1] = sszvec[k];
    }
    std::vector<int> rifvec(2 * mpisize, 0);
    SAFE_FUNC_EVAL( MPI_Alltoall( (void*)&(sifvec[0]), 2, MPI_INT, (void*)&(rifvec[0]), 2,
                      MPI_INT, MPI_COMM_WORLD ) );
    rszvec.resize(mpisize,0);
    for (int k = 0; k < mpisize; ++k) {
        _rnbvec[k] = rifvec[2 * k];
        rszvec[k] = rifvec[2 * k + 1];
    }
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::makeBufReqs(std::vector<int>& rszvec,
                                            std::vector<int>& sszvec) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::makeBufReqs");
#endif
    int mpisize = getMPISize();
    for (int k = 0; k < mpisize; ++k) {
        _rbufvec[k].resize(rszvec[k]);
    }
    for (int k = 0; k < mpisize; ++k) {
        // TODO (Austin): Is there a problem if rszvec or sszvec is 0?
        SAFE_FUNC_EVAL( MPI_Irecv((void *)&(_rbufvec[k][0]), rszvec[k], MPI_BYTE, k, 0,
                      MPI_COMM_WORLD, &_reqs[2 * k] ) );
        SAFE_FUNC_EVAL( MPI_Isend((void *)&(_sbufvec[k][0]), sszvec[k], MPI_BYTE, k, 0,
                      MPI_COMM_WORLD, &_reqs[2 * k + 1] ) );
        _kbytes_received += rszvec[k] / 1024;
        _kbytes_sent += sszvec[k] / 1024;
    }
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::strs2vec(std::vector<std::ostringstream *>& ossvec) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::strs2vec");
#endif
    int mpisize = getMPISize();
    for (int k = 0; k < mpisize; ++k) {
        std::string tmp( ossvec[k]->str() );
        _sbufvec[k].clear();
        _sbufvec[k].insert(_sbufvec[k].end(), tmp.begin(), tmp.end());
    }
    for (int k = 0; k < mpisize; ++k) {
        delete ossvec[k];
        ossvec[k] = NULL;
    }
    return 0;
}


//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::insert(Key key, Data& dat) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::insert");
#endif
    int mpirank = getMPIRank();
    CHECK_TRUE( _prtn.owner(key) == mpirank); //LEXING: VERY IMPORTANT
    _lclmap[key] = dat;
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
Data& ParVec<Key,Data,Partition>::access(Key key) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::access");
#endif
    typename std::map<Key,Data>::iterator mi = _lclmap.find(key);
    CHECK_TRUE(mi != _lclmap.end());
    return mi->second;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
std::pair<bool, Data&> ParVec<Key,Data,Partition>::contains(Key key) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::contains");
#endif
    typename std::map<Key,Data>::iterator mi = _lclmap.find(key);
    bool found = (mi  != _lclmap.end());
    if (found) {
      return std::pair<bool, Data&>(found, mi->second);
    }
    Data dummy;
    return std::pair<bool, Data&>(found, dummy);
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::getBegin(int (*e2ps)(Key, Data&, std::vector<int>&),
                                         const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::getBegin");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    resetVecs();
    _sbufvec.resize(mpisize);
    _rbufvec.resize(mpisize);
    _reqs = new MPI_Request[2 * mpisize];
    _stats = new MPI_Status[2 * mpisize];
    std::vector<std::ostringstream*> ossvec(mpisize);
    for (int k = 0; k < mpisize; ++k)        {
        ossvec[k] = new std::ostringstream();
    }

    // 1. serialize
    for (typename std::map<Key,Data>::iterator mi = _lclmap.begin();
        mi!=_lclmap.end(); ++mi) {
        Key key = mi->first;
        const Data& dat = mi->second;
        if (_prtn.owner(key) == mpirank) {
            std::vector<int> pids;
            int res = (*e2ps)(mi->first, mi->second, pids);
            for (int i = 0; i < pids.size(); ++i) {
                int k = pids[i];
                if (k != mpirank) { //DO NOT SEND TO MYSELF
                    SAFE_FUNC_EVAL( serialize(key, *(ossvec[k]), mask) );
                    SAFE_FUNC_EVAL( serialize(dat, *(ossvec[k]), mask) );
                    _snbvec[k]++; //LEXING: VERY IMPORTANT
                }
            }
        }
    }

    // to vector
    strs2vec(ossvec);

    // 2. all the sendsize of the message
    std::vector<int> sszvec;
    std::vector<int> rszvec;
    getSizes(rszvec, sszvec);

    // 3. allocate space, send and receive
    makeBufReqs(rszvec, sszvec);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::getBegin(std::vector<Key>& keyvec,
                                         const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::getBegin");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    resetVecs();
    _sbufvec.resize(mpisize);
    _rbufvec.resize(mpisize);
    _reqs = new MPI_Request[2 * mpisize];
    _stats = new MPI_Status[2 * mpisize];

    // 1. go through the keyvec to partition them among other procs
    std::vector< std::vector<Key> > skeyvec(mpisize);
    for (int i = 0; i < keyvec.size(); ++i) {
        Key key = keyvec[i];
        int owner = _prtn.owner(key);
        if (owner != mpirank) {
            skeyvec[owner].push_back(key);
        }
    }

    // 2. setdn receive size of keyvec
    std::vector<int> sszvec(mpisize);
    std::vector<int> rszvec(mpisize);
    for (int k = 0; k < mpisize; ++k) {
        sszvec[k] = skeyvec[k].size();
    }
    SAFE_FUNC_EVAL( MPI_Alltoall( (void*)&(sszvec[0]), 1, MPI_INT, (void*)&(rszvec[0]), 1,
                       MPI_INT, MPI_COMM_WORLD ) );

    // 3. allocate space for the keys, send and receive
    std::vector< std::vector<Key> > rkeyvec(mpisize);
    for (int k = 0; k < mpisize; ++k) {
        rkeyvec[k].resize(rszvec[k]);
    }

    MPI_Request *reqs = new MPI_Request[2 * mpisize];
    MPI_Status  *stats = new MPI_Status[2 * mpisize];
    for (int k = 0; k < mpisize; ++k) {
        SAFE_FUNC_EVAL( MPI_Irecv( (void*)&(rkeyvec[k][0]),
                                   rszvec[k] * sizeof(Key), MPI_BYTE,
                                   k, 0, MPI_COMM_WORLD, &reqs[2 * k] ) );
        SAFE_FUNC_EVAL( MPI_Isend( (void*)&(skeyvec[k][0]),
                                   sszvec[k] * sizeof(Key), MPI_BYTE,
                                   k, 0, MPI_COMM_WORLD, &reqs[2 * k + 1] ) );
    }
    SAFE_FUNC_EVAL( MPI_Waitall(2 * mpisize, &reqs[0], &stats[0]) );
    delete[] reqs;
    delete[] stats;

    skeyvec.clear(); // save space

    // 4. prepare the streams
    std::vector<std::ostringstream*> ossvec(mpisize);
    for (int k = 0; k < mpisize; ++k) {
        ossvec[k] = new std::ostringstream();
    }
    for (int k = 0; k < mpisize; ++k) {
        for (int g = 0; g < rkeyvec[k].size(); g++) {
            Key curkey = rkeyvec[k][g];
            typename std::map<Key, Data>::iterator mi = _lclmap.find(curkey);
            CHECK_TRUE( mi != _lclmap.end() );
            CHECK_TRUE( _prtn.owner(curkey) == mpirank );
            Key key = mi->first;
            const Data& dat = mi->second;
            SAFE_FUNC_EVAL( serialize(key, *(ossvec[k]), mask) );
            SAFE_FUNC_EVAL( serialize(dat, *(ossvec[k]), mask) );
            _snbvec[k]++; //LEXING: VERY IMPORTANT
        }
    }
    // to vector
    strs2vec(ossvec);

    // 5. all the sendsize of the message
    getSizes(rszvec, sszvec);

    // 6. allocate space, send and receive
    makeBufReqs(rszvec, sszvec);
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::getEnd( const std::vector<int>& mask ) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::getEnd");
#endif
    int mpisize = getMPISize();
    //LEXING: SEPARATE HERE
    SAFE_FUNC_EVAL( MPI_Waitall(2 * mpisize, &(_reqs[0]), &(_stats[0])) );
    delete[] _reqs;
    delete[] _stats;
    _sbufvec.clear(); //save space
    // Write back to stream.
    std::vector<std::istringstream*> issvec(mpisize);
    for (int k = 0; k < mpisize; ++k) {
        issvec[k] = new std::istringstream();
    }
    for (int k = 0; k < mpisize; ++k) {
        std::string tmp(_rbufvec[k].begin(), _rbufvec[k].end());
        issvec[k]->str(tmp);
    }
    _rbufvec.clear(); //save space
  
    for (int k = 0; k < mpisize; ++k) {
        for (int i = 0; i < _rnbvec[k]; ++i) {
            Key key;
            deserialize(key, *(issvec[k]), mask);
            typename std::map<Key, Data>::iterator mi = _lclmap.find(key);
            if (mi == _lclmap.end()) { //do not exist
                Data dat;
                deserialize(dat, *(issvec[k]), mask);
                _lclmap[key] = dat;
            } else { //exist already
                deserialize(mi->second, *(issvec[k]), mask);
            }
        }
    }
    for (int k = 0; k < mpisize; ++k) {
        delete issvec[k];
        issvec[k] = NULL;
    }
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::putBegin(std::vector<Key>& keyvec,
                                         const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::putBegin");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);
    resetVecs();
    _sbufvec.resize(mpisize);
    _rbufvec.resize(mpisize);
    _reqs = new MPI_Request[2 * mpisize];
    _stats = new MPI_Status[2 * mpisize];
    std::vector<std::ostringstream*> ossvec(mpisize);
    for (int k = 0; k < mpisize; ++k) {
        ossvec[k] = new std::ostringstream();
    }

    // 1. Go through the keyvec to partition them among other procs
    for (int i = 0; i < keyvec.size(); ++i) {
        Key key = keyvec[i];
        int k = _prtn.owner(key); //the owner
        if (k != mpirank) {
            typename std::map<Key, Data>::iterator mi = _lclmap.find(key);
            CHECK_TRUE( mi != _lclmap.end() );
            CHECK_TRUE( key == mi->first );
            Data& dat = mi->second;
            SAFE_FUNC_EVAL( serialize(key, *(ossvec[k]), mask) );
            SAFE_FUNC_EVAL( serialize(dat, *(ossvec[k]), mask) );
            _snbvec[k]++;
        }
    }

    // 2. to vector
    strs2vec(ossvec);

    // 3. get size
    std::vector<int> sszvec;
    std::vector<int> rszvec;
    getSizes(rszvec, sszvec);

    // 4. allocate space, send and receive
    makeBufReqs(rszvec, sszvec);

    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::putEnd( const std::vector<int>& mask ) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::putEnd");
#endif
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    //LEXING: SEPARATE HERE
    SAFE_FUNC_EVAL( MPI_Waitall(2*mpisize, &(_reqs[0]), &(_stats[0])) );
    delete[] _reqs;
    delete[] _stats;
    _sbufvec.clear(); //save space
    //5. go thrw the messages and write back
    std::vector<std::istringstream*> issvec(mpisize);
    for (int k = 0; k < mpisize; ++k) {
        issvec[k] = new std::istringstream();
    }
    for (int k = 0; k < mpisize; ++k) {
        std::string tmp(_rbufvec[k].begin(), _rbufvec[k].end());
        issvec[k]->str(tmp);
    }
    _rbufvec.clear(); //save space
    for (int k = 0; k < mpisize; ++k) {
        for (int i = 0; i < _rnbvec[k]; ++i) {
            Key key;
            deserialize(key, *(issvec[k]), mask);
            CHECK_TRUE( _prtn.owner(key) == mpirank );
            typename std::map<Key,Data>::iterator mi = _lclmap.find(key);
            CHECK_TRUE( mi!=_lclmap.end() );
            deserialize(mi->second, *(issvec[k]), mask);
        }
    }
    for (int k = 0; k < mpisize; ++k) {
        delete issvec[k];
        issvec[k] = NULL;
    }
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::expand(std::vector<Key>& keyvec) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::expand");
#endif
    Data dummy;
    for (int i = 0; i<keyvec.size(); ++i) {
        Key key = keyvec[i];
        typename std::map<Key, Data>::iterator mi = _lclmap.find(key);
        if (mi == _lclmap.end()) {
            _lclmap[key] = dummy;
        }
    }
    return 0;
}

//--------------------------------------------
template <class Key, class Data, class Partition>
int ParVec<Key,Data,Partition>::discard(std::vector<Key>& keyvec) {
#ifndef RELEASE
    CallStackEntry entry("ParVec::discard");
#endif
    int mpirank = getMPIRank();
    for (int i = 0; i < keyvec.size(); ++i) {
        Key key = keyvec[i];
        if (_prtn.owner(key) != mpirank) {
            _lclmap.erase(key);
        }
    }
    return 0;
}

//-------------------
template<class Key, class Data, class Partition>
int serialize(const ParVec<Key,Data,Partition>& pv, std::ostream& os,
              const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("serialize");
#endif
    serialize(pv._lclmap, os, mask);
    serialize(pv._prtn, os, mask);
    return 0;
}

template<class Key, class Data, class Partition>
int deserialize(ParVec<Key,Data,Partition>& pv, std::istream& is,
                const std::vector<int>& mask) {
#ifndef RELEASE
    CallStackEntry entry("deserialize");
#endif
    deserialize(pv._lclmap, is, mask);
    deserialize(pv._prtn, is, mask);
    return 0;
}

#endif
