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
#include "parallel.hpp"

#include "commoninc.hpp"

#include <iostream>
#include <sstream>
#include <string>

#define MAX_FILE_NAME_LENGTH 128

//---------------------------------------------------------
int SeparateRead(std::string name, std::istringstream& is) {
#ifndef RELEASE
    CallStackEntry entry("Separate_Read");
#endif
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    char filename[MAX_FILE_NAME_LENGTH];
    sprintf(filename, "data/%s_%d_%d", name.c_str(), mpirank, mpisize);
    std::ifstream fin(filename);
    if (fin.fail()) {
	std::cout << "WARNING: failed to open input file stream: " << filename << std::endl;
	return 0;
    }
    is.str( std::string(std::istreambuf_iterator<char>(fin),
			std::istreambuf_iterator<char>()) );
    fin.close();
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//---------------------------------------------------------
int SeparateWrite(std::string name, std::ostringstream& os) {
#ifndef RELEASE
    CallStackEntry entry("Separate_Write");
#endif
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    char filename[MAX_FILE_NAME_LENGTH];
    sprintf(filename, "data/%s_%d_%d", name.c_str(), mpirank, mpisize);
    std::ofstream fout(filename);
    if (fout.fail()) {
	std::cerr << "failed to open output file stream: " << filename << std::endl;
	return -1;
    }
    fout<<os.str();
    fout.close();
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//---------------------------------------------------------
int SharedRead(std::string name, std::istringstream& is) {
#ifndef RELEASE
    CallStackEntry entry("Shared_Read");
#endif
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    std::vector<char> tmpstr;
    if (mpirank == 0) {
	char filename[MAX_FILE_NAME_LENGTH];
	sprintf(filename, "data/%s", name.c_str());
	std::cerr << filename << std::endl;
	std::ifstream fin(filename);
	if (fin.fail()) {
	    fprintf(stderr, "failed to open input file stream (%s)\n", filename);
	    return -1;
	}

	tmpstr.insert(tmpstr.end(), std::istreambuf_iterator<char>(fin),
		      std::istreambuf_iterator<char>());
	fin.close();
	int size = tmpstr.size();
	SAFE_FUNC_EVAL( MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD) );
	SAFE_FUNC_EVAL( MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD) );
    } else {
	int size;
	SAFE_FUNC_EVAL( MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD) );
	tmpstr.resize(size);
	SAFE_FUNC_EVAL( MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD) );
    }
    is.str( std::string(tmpstr.begin(), tmpstr.end()) );
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}

//---------------------------------------------------------
int SharedWrite(std::string name, std::ostringstream& os) {
#ifndef RELEASE
    CallStackEntry entry("Shared_Write");
#endif
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    int mpirank, mpisize;
    getMPIInfo(&mpirank, &mpisize);

    //
    if (mpirank == 0) {
	char filename[MAX_FILE_NAME_LENGTH];
	sprintf(filename, "data/%s", name.c_str());
	std::cerr << filename << std::endl;
	std::ofstream fout(filename);
	if (fout.fail()) {
	    fprintf(stderr, "failed to open output file stream (%s)\n", filename);
	    return -1;
	}
	fout << os.str();
	fout.close();
    }
    SAFE_FUNC_EVAL( MPI_Barrier(MPI_COMM_WORLD) );
    return 0;
}
