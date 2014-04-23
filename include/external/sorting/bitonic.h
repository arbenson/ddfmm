#ifndef __BITONIC_H_
#define __BITONIC_H_

#include <cassert>
#include <vector>
#include "mpi.h"

// #include "Quad.h"
#include "dtypes.h"

//#include "colors.h"


#define KEEP_HIGH 100
#define KEEP_LOW  101


bool isPowerOfTwo(unsigned int n) {
  return !(n & (n - 1)) && n;
}

// compute the next highest power of 2 of 32-bit v
int getNextHighestPowerOfTwo(unsigned int n) {
  unsigned int v = n;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

// compute the prev highest power of 2 of 32-bit v
int getPrevHighestPowerOfTwo(unsigned int n) {
  unsigned int v = n;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return  v >> 1;
}

/** 
 * @brief Splits a communication group into two, the first having a power of 2
 * number of processors and the other having the remainder.
 * 
 * @param orig_comm    The comm group that needs to be split.
 * @param new_comm     The new comm group.
 */
unsigned int splitCommBinary( MPI_Comm orig_comm, MPI_Comm *new_comm) {
  int npes, rank;
  
  MPI_Group  orig_group, new_group;

  MPI_Comm_size(orig_comm, &npes);
  MPI_Comm_rank(orig_comm, &rank);

  unsigned int splitterRank = getPrevHighestPowerOfTwo(npes);

  int *ranksAsc, *ranksDesc;
  //Determine sizes for the 2 groups 
  ranksAsc = new int[splitterRank];
  ranksDesc = new int[( npes - splitterRank)];

  int numAsc = 0;
  int numDesc = ( npes - splitterRank - 1);

  //This is the main mapping between old ranks and new ranks.
  for(int i=0; i<npes; i++) {
    if(i < splitterRank) {
      ranksAsc[numAsc] = i;
        numAsc++;
    }else {
      ranksDesc[numDesc] = i;
        numDesc--;
    }
  }//end for i

  MPI_Comm_group(orig_comm, &orig_group);

  /* Divide tasks into two distinct groups based upon rank */
  if (rank < splitterRank) {
    MPI_Group_incl(orig_group, splitterRank, ranksAsc, &new_group);
  }else {
    MPI_Group_incl(orig_group, (npes-splitterRank), ranksDesc, &new_group);
  }
  
  MPI_Comm_create(orig_comm, new_group, new_comm);

  delete [] ranksAsc;
  delete [] ranksDesc;

  return splitterRank;
}

unsigned int splitCommBinaryNoFlip( MPI_Comm orig_comm, MPI_Comm *new_comm) {
  int npes, rank;
  
  MPI_Group  orig_group, new_group;

  MPI_Comm_size(orig_comm, &npes);
  MPI_Comm_rank(orig_comm, &rank);

  unsigned int splitterRank = getPrevHighestPowerOfTwo(npes);

  int *ranksAsc, *ranksDesc;
  //Determine sizes for the 2 groups 
  ranksAsc = new int[splitterRank];
  ranksDesc = new int[( npes - splitterRank)];

  int numAsc = 0;
  int numDesc = 0; //( npes - splitterRank - 1);

  //This is the main mapping between old ranks and new ranks.
  for(int i=0; i<npes; i++) {
    if(i < splitterRank) {
      ranksAsc[numAsc] = i;
        numAsc++;
    }else {
      ranksDesc[numDesc] = i;
        numDesc++;
    }
  }//end for i

  MPI_Comm_group(orig_comm, &orig_group);

  /* Divide tasks into two distinct groups based upon rank */
  if (rank < splitterRank) {
    MPI_Group_incl(orig_group, splitterRank, ranksAsc, &new_group);
  }else {
    MPI_Group_incl(orig_group, (npes-splitterRank), ranksDesc, &new_group);
  }
  
  MPI_Comm_create(orig_comm, new_group, new_comm);

  delete [] ranksAsc;
  delete [] ranksDesc;

  return splitterRank;
}


/** 
 * @brief Merges lists A, and B, retaining either the low or the High in list A.      
 * 
 * @param listA   	Input list, and where the output is stored.
 * @param listB   	Second input list.
 * @param KEEP_WHAT 	determines whether to retain the High or the low values
 * 			from A and B. One of KEEP_HIGH or KEEP_LOW.
 *
 * Merging the two lists when their sizes are not the same is a bit involved.
 * The major condition that needs to be used is that all elements that are less
 * than max(min(A), min(B)) are retained by the KEEP_LOW processor, and
 * similarly all elements that are larger larger than min(max(A), max(B)) are
 * retained by the KEEP_HIGH processor.
 *
 * The reason for this is that, on the Keep_Low side,
 *
 *   max(min(A), min(B)) > min(A) > max(A-)
 *
 * and similarly on the Keep_high side,
 *  
 *   min(max(A), max(B)) < max(A) < min(A+)
 *
 * which guarantees that the merged lists remain bitonic.
 */

template <typename T>
void MergeLists( std::vector<T> &listA, std::vector<T> &listB, int KEEP_WHAT) {
  unsigned int list_size;
  std::vector<T> scratch_list;
  
  T _low, _high;

  assert(listA.size());
  assert(listB.size());
  
  _low  = (listA[0]>listB[0])?listA[0]:listB[0];
  _high = (listA[listA.size()-1]<listB[listB.size()-1])?listA[listA.size()-1]:listB[listB.size()-1];

  // We will do a full merge first ... 
  list_size = (listA.size() + listB.size());
  scratch_list.resize(list_size);
  
  unsigned int  index1 = 0;
  unsigned int  index2 = 0;  

  for (unsigned int i = 0; i < list_size; i++) {
    if ( (index1 < listA.size()) && ( (listA[index1] <= listB[index2]) || !(index2 < listB.size()) ) ) {
      scratch_list[i] = listA[index1++];
    } else {
      scratch_list[i] = listB[index2++];
    }
  }

  listA.clear();
  if ( KEEP_WHAT == KEEP_LOW ) {
    int ii=0;
    while ( ((scratch_list[ii] < _low) || (ii < list_size/2) ) && (scratch_list[ii] <= _high) )
      listA.push_back(scratch_list[ii++]);

  } else {
    int ii=list_size-1;

    while ( (scratch_list[ii] >= _low) && ( (ii >= list_size/2) || (scratch_list[ii] > _high) ) )
      listA.insert(listA.begin(), scratch_list[ii--]);
  }

  scratch_list.clear();
}

/********************************************************************/
/*
 * which_keys is one of KEEP_HIGH or KEEP_LOW
 * partner    is the processor with which to Merge and Split.
 *
 */

template <typename T>
void MergeSplit( std::vector<T> &local_list, int which_keys, int partner, MPI_Comm  comm) {

  MPI_Status status;
  unsigned int send_size = local_list.size();
  unsigned int recv_size = local_list.size();

  // first communicate how many you will send and how many you will receive ...

  MPI_Sendrecv( &send_size , 1, MPI_UNSIGNED, partner, 0, &recv_size, 1, MPI_UNSIGNED, partner, 0, comm, &status);

  std::vector<T> temp_list( recv_size );

  MPI_Sendrecv( &(*local_list.begin()) , send_size, par::Mpi_datatype<T>::value(), partner, 1, &(*temp_list.begin()), recv_size, par::Mpi_datatype<T>::value(), partner, 1, comm, &status);

  MergeLists<T>(local_list, temp_list, which_keys);

  temp_list.clear();
} /* Merge_split */




/********************************************************************/
template <typename T>
void Par_bitonic_sort_incr( std::vector<T> &local_list, int proc_set_size, MPI_Comm  comm ) {
  unsigned  eor_bit;
  int       proc_set_dim;
  int       stage;
  int       partner;
  int       my_rank, npes;

  MPI_Comm_rank(comm, &my_rank);

  proc_set_dim = 0;
  int x = proc_set_size;
  while (x > 1) {
    x = x >> 1;
    proc_set_dim++;
  }

  eor_bit = 1 << (proc_set_dim - 1);
  for (stage = 0; stage < proc_set_dim; stage++) {
    partner = my_rank ^ eor_bit;
    
    if (my_rank < partner)
      MergeSplit<T> ( local_list,  KEEP_LOW, partner, comm);
    else
      MergeSplit<T> ( local_list, KEEP_HIGH, partner, comm);
    
    eor_bit = eor_bit >> 1;
  }
}  /* Par_bitonic_sort_incr */


/********************************************************************/
template <typename T>
void Par_bitonic_sort_decr( std::vector<T> &local_list, int proc_set_size, MPI_Comm  comm) {
  unsigned  eor_bit;
  int       proc_set_dim;
  int       stage;
  int       partner;
  int       my_rank, npes;

  MPI_Comm_rank(comm, &my_rank);
  
  proc_set_dim = 0;
  int x = proc_set_size;
  while (x > 1) {
    x = x >> 1;
    proc_set_dim++;
  }

  eor_bit = 1 << (proc_set_dim - 1);
  for (stage = 0; stage < proc_set_dim; stage++) {
    partner = my_rank ^ eor_bit;
    
    if (my_rank > partner)
      MergeSplit<T> ( local_list,  KEEP_LOW, partner, comm);
    else
      MergeSplit<T> ( local_list, KEEP_HIGH, partner, comm);

    eor_bit = eor_bit >> 1;
  }

} /* Par_bitonic_sort_decr */

template <typename T>
void Par_bitonic_merge_incr( std::vector<T> &local_list, MPI_Comm  comm ) {
  int       partner;
  int       rank, npes;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &npes);

  unsigned int num_left  = getPrevHighestPowerOfTwo(npes);
  unsigned int num_right = npes - num_left;
  
  // 1, Do merge between the k right procs and the highest k left procs.
  if ( (rank < num_left) && ( rank >= (num_left - num_right) ) ) {
    partner = rank + num_right;
    // std::cout << RED << rank << YLW" --> "GRN << partner << NRM << std::endl;
    MergeSplit<T> ( local_list,  KEEP_LOW, partner, comm);
  } else if (rank >= num_left) {
    partner = rank - num_right;
    // std::cout << BLU << rank << CYN" <-- "MAG << partner << NRM << std::endl;
    MergeSplit<T> ( local_list,  KEEP_HIGH, partner, comm);
  }
}


template <typename T>
void bitonicSort_binary(std::vector<T> & in, MPI_Comm comm) {
  int       	    proc_set_size;
  unsigned int	    and_bit;
  int       	rank;
  int       	npes;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &rank);

  for (proc_set_size = 2, and_bit = 2;
      proc_set_size <= npes;
      proc_set_size = proc_set_size*2, 
      and_bit = and_bit << 1) {

    if ((rank & and_bit) == 0)
      Par_bitonic_sort_incr<T>( in, proc_set_size, comm);
    else
      Par_bitonic_sort_decr<T>( in, proc_set_size, comm);
    MPI_Barrier(comm);
  }
}

template <typename T>
void bitonicSort(std::vector<T> & in, MPI_Comm comm) {
  int       	rank;
  int       	npes;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &rank);

  // local sort first ...
  std::sort (in.begin(), in.end());

  // check if npes is a power of two ...
  bool isPower = !(npes & npes - 1) && npes ;

  if ( isPower ) {
    bitonicSort_binary<T>(in, comm);
  } else {
    MPI_Comm new_comm;

    // Since npes is not a power of two, we shall split the problem in two ...
    //
    // 1. Create 2 comm groups ... one for the 2^d portion and one for the
    // remainder.
    unsigned int splitter = splitCommBinary(comm, &new_comm);
    
    if ( rank < splitter)
      bitonicSort_binary<T>(in, new_comm);
    else
      bitonicSort<T>(in, new_comm);
 
    // 3. Do a special merge of the two segments. (original comm).
    Par_bitonic_merge_incr( in, comm );
   
    splitter = splitCommBinaryNoFlip(comm, &new_comm);
    
    // 4. Now a final sort on the segments.
    if ( rank < splitter)
      bitonicSort_binary<T>(in, new_comm);
    else
      bitonicSort<T>(in, new_comm);
  }
}


#endif
