#include "parallel.hpp"

using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;
using std::endl;

#define MAX_FILE_NAME_LENGTH 128

//---------------------------------------------------------
int Separate_Read(string name, istringstream& is)
{
#ifndef RELEASE
    CallStackEntry entry("Separate_Read");
#endif
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    char filename[MAX_FILE_NAME_LENGTH];
    sprintf(filename, "data/%s_%d_%d", name.c_str(), mpirank, mpisize);  
    cerr << filename << endl;
    ifstream fin(filename);
    if (fin.fail()) {
	fprintf(stderr, "failed to open input file stream (%s)\n", filename);
    }
    // TODO (Austin): We should probably exit in this case
    is.str( string(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>()) );
    fin.close();
    //
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//---------------------------------------------------------
int Separate_Write(string name, ostringstream& os)
{
#ifndef RELEASE
    CallStackEntry entry("Separate_Write");
#endif
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    //
    char filename[MAX_FILE_NAME_LENGTH];
    sprintf(filename, "data/%s_%d_%d", name.c_str(), mpirank, mpisize);
    cerr << filename << endl;
    ofstream fout(filename);
    if (fout.fail()) {
	fprintf(stderr, "failed to open output file stream (%s)\n", filename);
    }
    // TODO (Austin): We should probably exit in this case
    fout<<os.str();
    fout.close();
    //
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//---------------------------------------------------------
int Shared_Read(string name, istringstream& is)
{
#ifndef RELEASE
    CallStackEntry entry("Shared_Read");
#endif
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    vector<char> tmpstr;
    if(mpirank == 0) {
	char filename[MAX_FILE_NAME_LENGTH];
	sprintf(filename, "data/%s", name.c_str());
	cerr << filename << endl;
	ifstream fin(filename);
	if (fin.fail()) {
	    // TODO (Austin): We should probably exit in this case
	    fprintf(stderr, "failed to open input file stream (%s)\n", filename);
	}

	tmpstr.insert(tmpstr.end(), std::istreambuf_iterator<char>(fin),
		      std::istreambuf_iterator<char>());
	fin.close();
	int size = tmpstr.size();
	iC( MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD) );
	iC( MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD) );
    } else {
	int size;
	iC( MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD) );
	tmpstr.resize(size);
	iC( MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD) );
    }
    is.str( string(tmpstr.begin(), tmpstr.end()) );
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//---------------------------------------------------------
int Shared_Write(string name, ostringstream& os)
{
#ifndef RELEASE
    CallStackEntry entry("Shared_Write");
#endif
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    //
    if(mpirank == 0) {
	char filename[MAX_FILE_NAME_LENGTH];
	sprintf(filename, "data/%s", name.c_str());
	cerr << filename << endl;
	ofstream fout(filename);
	if (fout.fail()) {
	    // TODO (Austin): We should probably exit in this case
	    fprintf(stderr, "failed to open output file stream (%s)\n", filename);
	}
	fout << os.str();
	fout.close();
    }
    iC( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}
