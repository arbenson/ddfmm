Distributed Directional Fast Multipole Method
=====
### Austin R. Benson, Jack Poulson, Lexing Ying

This code provides a distributed memory implementation of the directional Fast Multipole Method.
Given a set of N points p_1, ..., p_N and N densities f_1, ..., f_N, we want to compute the potentials
u_1, ..., u_N defined by

u_i = G(p_i, p_1)f_1 + ... + G(p_i, p_N)f_N,

where G(p_i, p_j) = exp(i|p_i - p_j|) / |p_i - p_j| is the Green's function of the Helmholtz equation.
Our code also supports the kernel G(p_i, p_j) = exp(i|p_i - p_j|).

The code is liencensed under GPLv3.  Please see the COPYING file for details of the license.

Pre-computation and data generation
-----
The files in the matlab directory are tools for the pre-computation and data generation.
The pre-computation contains all of the translation matrices for the low and high frequency regimes.
The generated data are the points and densities sampled from one of the geometries.
We provide sphere, F16, and submarine geometries.

# Generate the translation matrices for high and low frequency
cd matlab
matlab < aug3d_script.m

# Generate the points, densities, and partitions for different geometries
cd matlab
matlab < data_script.m


Building and running the code
-----
The code requires the following libraries:
* mpi
* fftw
* blas
* lapack

There is also support for MKL.  You need to compile with the variable MKL defined.
See the file makeinc/lonestar for an example.

To build the program:
1. Define the environment variable HOST and create the file corresponding file makeinc/${HOST}.
2. make tt

Running the code
-----
Edit and use the file run.sh.

Contact
--------
For questions, suggestions, and bug reports, please email Austin Benson: arbenson AT stanford DOT edu.
