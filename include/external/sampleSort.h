#ifndef _SAMPLE_SORT_H_
#define _SAMPLE_SORT_H_ 

#include "bitonic.h"

template <typename T>
int UpperBound (unsigned int p,const T * splitt,unsigned int itr, const T & elem)
{
  if (itr >= p) {
    return p;
  }
  while (itr < p){
    if (elem <= splitt[itr]) {
      return itr;
    } else {
      itr = itr + 1;
    }
  }//end while
  return itr;
}//end function

int splitCommUsingSplittingRank(int splittingRank, MPI_Comm* new_comm,
    MPI_Comm comm) {

  MPI_Group  orig_group, new_group;
  int size;
  int rank;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int* ranksActive = new int[splittingRank];
  int* ranksIdle = new int[size - splittingRank];

  for(int i = 0; i < splittingRank; i++) {
    ranksActive[i] = i;
  }

  for(int i = splittingRank; i < size; i++) {
    ranksIdle[i - splittingRank] = i;
  }

  /* Extract the original group handle */
  MPI_Comm_group(comm, &orig_group);

  /* Divide tasks into two distinct groups based upon rank */
  if (rank < splittingRank) {
    MPI_Group_incl(orig_group, splittingRank, ranksActive, &new_group);
  }else {
    MPI_Group_incl(orig_group, (size - splittingRank), ranksIdle, &new_group);
  }

  /* Create new communicator */
  MPI_Comm_create(comm, new_group, new_comm);

  delete [] ranksActive;
  ranksActive = NULL;

  delete [] ranksIdle;
  ranksIdle = NULL;

  return 1;
}//end function

template<typename T>
unsigned int defaultWeight(const T *a){
  return 1;
}

template<typename T>
int partitionW(std::vector<T>& nodeList,  MPI_Comm comm){

  int npes;

  MPI_Comm_size(comm, &npes);

  int rank;

  MPI_Comm_rank(comm, &rank);

  MPI_Request request;
  MPI_Status status;
  const bool nEmpty = nodeList.empty();

  int  off1= 0, off2= 0, localWt= 0, totalWt = 0;

  int* wts = NULL;
  int* lscn = NULL;
  int nlSize = nodeList.size();
  if(nlSize) {
    wts = new int[nlSize];
    assert(wts);

    lscn= new int[nlSize]; 
    assert(lscn);
  }

  // First construct arrays of id and wts.
  for (int i = 0; i < nlSize; i++) {
    wts[i] = 1; // (*getWeight)( &(nodeList[i]) );
    localWt += wts[i];
  }

  // compute the total weight of the problem ...
  MPI_Allreduce (&localWt, &totalWt, 1, MPI_INT, MPI_SUM, comm);

  // perform a local scan on the weights first ...
  int zero = 0;
  if(!nEmpty) {
    lscn[0]=wts[0];
    for (int i = 1; i < nlSize; i++) {
      lscn[i] = wts[i] + lscn[i-1];
    }//end for
    // now scan with the final members of 
    MPI_Scan (lscn+nlSize-1, &off1, 1, MPI_INT, MPI_SUM, comm ); 
  } else{
    MPI_Scan (&zero, &off1, 1, MPI_INT, MPI_SUM, comm ); 
  }

  // communicate the offsets ...
  if (rank < (npes-1)){
    MPI_Issend ( &off1, 1, MPI_INT, rank+1, 0, comm, &request );
  }
  if (rank){
    MPI_Recv ( &off2, 1, MPI_INT, rank-1, 0, comm, &status );
  }
  else{
    off2 = 0; 
  }

  // add offset to local array
  for (int i = 0; i < nlSize; i++) {
    lscn[i] = lscn[i] + off2;       // This has the global scan results now ...
  }//end for

  int * sendSz = new int [npes];
  assert(sendSz);

  int * recvSz = new int [npes];
  assert(recvSz);

  int * sendOff = new int [npes]; 
  assert(sendOff);
  sendOff[0] = 0;

  int * recvOff = new int [npes]; 
  assert(recvOff);
  recvOff[0] = 0;

  // compute the partition offsets and sizes so that All2Allv can be performed.
  // initialize ...

  for (int i = 0; i < npes; i++) {
    sendSz[i] = 0;
  }

  // Now determine the average load ...
  int npesLong = npes;
  int avgLoad = (totalWt/npesLong);

  int extra = (totalWt%npesLong);

  //The Heart of the algorithm....
  if(avgLoad > 0) {
    for (int i = 0; i < nlSize; i++) {
      if(lscn[i] == 0) {		
        sendSz[0]++;
      }else {
        int ind=0;
        if ( lscn[i] <= (extra*(avgLoad + 1)) ) {
          ind = ((lscn[i] - 1)/(avgLoad + 1));
        }else {
          ind = ((lscn[i] - (1 + extra))/avgLoad);
        }
        assert(ind < npes);
        sendSz[ind]++;
      }//end if-else
    }//end for 
  }else {
    sendSz[0]+= nlSize;
  }//end if-else

  if(rank < (npes-1)) {
    MPI_Status statusWait;
    MPI_Wait(&request, &statusWait);
  }

  // communicate with other procs how many you shall be sending and get how
  // many to recieve from whom.
  MPI_Alltoall (sendSz, 1, MPI_INT, recvSz, 1, MPI_INT, comm);

  int nn=0; // new value of nlSize, ie the local nodes.
  for (int i = 0; i < npes; i++) {
    nn += recvSz[i];
  }

  // compute offsets ...
  for (int i = 1; i < npes; i++) {
    sendOff[i] = sendOff[i-1] + sendSz[i-1];
    recvOff[i] = recvOff[i-1] + recvSz[i-1];
  }

  // allocate memory for the new arrays ...
  std::vector<T > newNodes(nn);

  // perform All2All  ... 
  T* nodeListPtr = NULL;
  T* newNodesPtr = NULL;
  if(!nodeList.empty()) {
    nodeListPtr = &(*(nodeList.begin()));
  }
  if(!newNodes.empty()) {
    newNodesPtr = &(*(newNodes.begin()));
  }
  MPI_Alltoallv (nodeListPtr, sendSz, sendOff, par::Mpi_datatype<T>::value(),
                            newNodesPtr, recvSz, recvOff, par::Mpi_datatype<T>::value(), comm);

  // reset the pointer ...
  nodeList = newNodes;
  newNodes.clear();

  // clean up...
  if(!nEmpty) {
    delete [] lscn;
    delete [] wts;
  }
  delete [] sendSz;
  sendSz = NULL;

  delete [] sendOff;
  sendOff = NULL;

  delete [] recvSz;
  recvSz = NULL;

  delete [] recvOff;
  recvOff = NULL;

  return 1;
}//end function

template<typename T>
int sampleSort(std::vector<T>& arr, std::vector<T> & SortedElem, MPI_Comm comm){ 

  int npes;

  MPI_Comm_size(comm, &npes);

  if (npes == 1) {
    std::cout <<" have to use seq. sort"
      <<" since npes = 1 . inpSize: "<<(arr.size()) <<std::endl;
    std::sort(arr.begin(), arr.end());
    SortedElem  = arr;
    return 1;
  } 

  std::vector<T>  splitters;
  std::vector<T>  allsplitters;

  int myrank;
  MPI_Comm_rank(comm, &myrank);

  int nelem = arr.size();
  int nelemCopy = nelem;
  int totSize;
  MPI_Allreduce(&nelemCopy, &totSize, 1, MPI_INT, MPI_SUM, comm);

  int npesLong = npes;
  const int FIVE = 5;

  if(totSize < (FIVE*npesLong)) {
    if(!myrank) {
      std::cout <<" Using bitonic sort since totSize < (5*(npes)). totSize: "
        <<totSize<<" npes: "<<npes <<std::endl;
    }
    partitionW<T>(arr, comm);

    SortedElem = arr; 
    MPI_Comm new_comm;
    if(totSize < npesLong) {
      if(!myrank) {
        std::cout<<" Input to sort is small. splittingComm: "
          <<npes<<" -> "<< totSize<<std::endl;
      }
      splitCommUsingSplittingRank(static_cast<int>(totSize), &new_comm, comm);
    } else {
      new_comm = comm;
    }

    if(!SortedElem.empty()) {
      bitonicSort<T>(SortedElem, new_comm);
    }

  }// end if

  //Re-part arr so that each proc. has atleast p elements.
  partitionW<T>(arr, comm);

  nelem = arr.size();

  std::sort(arr.begin(),arr.end());

  std::vector<T> sendSplits(npes-1);
  splitters.resize(npes);

  for(int i = 1; i < npes; i++)	 {
    sendSplits[i-1] = arr[i*nelem/npes];	
  }//end for i

  // sort sendSplits using bitonic ...
  bitonicSort<T>(sendSplits,comm);

  // All gather with last element of splitters.
  T* sendSplitsPtr = NULL;
  T* splittersPtr = NULL;
  if(sendSplits.size() > static_cast<unsigned int>(npes-2)) {
    sendSplitsPtr = &(*(sendSplits.begin() + (npes -2)));
  }
  if(!splitters.empty()) {
    splittersPtr = &(*(splitters.begin()));
  }
  MPI_Allgather (sendSplitsPtr, 1, par::Mpi_datatype<T>::value(), splittersPtr, 1, par::Mpi_datatype<T>::value(), comm);

  sendSplits.clear();

  int *sendcnts = new int[npes];
  assert(sendcnts);

  int * recvcnts = new int[npes];
  assert(recvcnts);

  int * sdispls = new int[npes];
  assert(sdispls);

  int * rdispls = new int[npes];
  assert(rdispls);

  for(int k = 0; k < npes; k++){
    sendcnts[k] = 0;
  }

  int k = 0;

  for (int j = 0; j < nelem; j++) {
    if (arr[j] <= splitters[k]) {
      sendcnts[k]++;
    } else{
      k = UpperBound<T>(npes-1, splittersPtr, k+1, arr[j]);
      if (k == (npes-1) ){
        //could not find any splitter >= arr[j]
        sendcnts[k] = (nelem - j);
        break;
      } else {
        assert(k < (npes-1));
        assert(splitters[k] >= arr[j]);
        sendcnts[k]++;
      }
    }//end if-else
  }//end for j

  MPI_Alltoall (sendcnts, 1, MPI_INT, recvcnts, 1, MPI_INT, comm);

  sdispls[0] = 0; rdispls[0] = 0;
  for (int j = 1; j < npes; j++){
    sdispls[j] = sdispls[j-1] + sendcnts[j-1];
    rdispls[j] = rdispls[j-1] + recvcnts[j-1];
  }

  int nsorted = rdispls[npes-1] + recvcnts[npes-1];
  SortedElem.resize(nsorted);

  T* arrPtr = NULL;
  T* SortedElemPtr = NULL;
  if(!arr.empty()) {
    arrPtr = &(*(arr.begin()));
  }
  if(!SortedElem.empty()) {
    SortedElemPtr = &(*(SortedElem.begin()));
  }
  MPI_Alltoallv   (arrPtr, sendcnts, sdispls, par::Mpi_datatype<T>::value(), 
      SortedElemPtr, recvcnts, rdispls, par::Mpi_datatype<T>::value(), comm);

  arr.clear();

  delete [] sendcnts;
  sendcnts = NULL;

  delete [] recvcnts;
  recvcnts = NULL;

  delete [] sdispls;
  sdispls = NULL;

  delete [] rdispls;
  rdispls = NULL;

  sort(SortedElem.begin(), SortedElem.end());

  return 1;
}//end function

#endif 
