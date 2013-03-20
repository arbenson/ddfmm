#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;

#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[]) {
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  printf("Hello, world! I am %d of %d\n", rank,size);
  printf("%d\n", argc);
  
  MPI_Finalize();
  return 0;
}

/*
int main(int argc, char *argv[]) {
  ifstream fin;
  fin.open("aa");
  //string str((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
  //cerr<<str<<endl;
  //ofstream fout("bb");
  //fout<<str;
  istringstream iss;
  iss.str( string((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>()) );
  cerr<<iss.str();
  return 0;
}
*/

/*
  int main ()
{
  ostringstream oss (ostringstream::out);
  oss << "This is a test\n";
  cout << oss.str();
  vector< ostringstream > gg;
  gg.push_back( ostringstream() );
  gg.push_back( ostringstream() );
  return 0;
}
*/
