#include "parallel.hpp"

using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;

//---------------------------------------------------------
int Separate_Read(string name, istringstream& is)
{
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);  //cerr<<filename<<endl;
  ifstream fin(filename);
  is.str( string(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>()) );
  fin.close();
  //
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  return 0;
}

//---------------------------------------------------------
int Separate_Write(string name, ostringstream& os)
{
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);
  ofstream fout(filename);
  fout<<os.str();
  fout.close();
  //
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  return 0;
}

//---------------------------------------------------------
int Shared_Read(string name, istringstream& is)
{
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  vector<char> tmpstr;
  if(mpirank==0) {
    ifstream fin(name.c_str());
    //string str(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
    //tmpstr.insert(tmpstr.end(), str.begin(), str.end());
    tmpstr.insert(tmpstr.end(), std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
    fin.close();
    int size = tmpstr.size();	//cerr<<size<<endl;
    iC( MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD) );
    iC( MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD) );
  } else {
    int size;
    iC( MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD) );
    tmpstr.resize(size);
    iC( MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD) );
  }
  is.str( string(tmpstr.begin(), tmpstr.end()) );
  //
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  return 0;
}

//---------------------------------------------------------
int Shared_Write(string name, ostringstream& os)
{
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  if(mpirank==0) {
    ofstream fout(name.c_str());
    fout<<os.str();
    fout.close();
  }
  iC( MPI_Barrier(MPI_COMM_WORLD) );
  return 0;
}
