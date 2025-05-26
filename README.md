# A PARALLEL JACOBI ITERATION SOLVER FOR THE POISSON PROBLEM
### S. Galati
### A.A. 2024-2025

## Code organization
The code is organized in three directories: `src`, `include`, `lib` and `test`. The first two directories contain the implementation of the laplace solver, organized as a small library, the lib folder contains the muParser dynamic library, while the test folder contains the Makefile, the datafile and a bashscript to make a scalability test.

Codes taken from the PACS course:
- `chrono.hpp`: utility to test execution time
- `mpi_utils.hpp`: utility to allor generality using MPI types

### Source code
Three main classes are used:
`Params`:
A struct holding the solver parameters.
- It implements a function for reading the parameters from a GetPot datafile, including some options passed at runtime.
- Notice that, in order to Broadcast efficiently minimizing communication, the strings are stored apart from the other parameters so that they can be sent with a single Broadcast after concatenating them.

`solMatrix`: A class to store a (local) parallel row-major dense matrix.
- Each localMatrix owns a certain number of rows (M_nrowsInternal) and performs local computations only on those lines. Moreover it stores 2 ghost matrices, one for bottom and one for top. For internal ranks, those rows will be used to receive the current solution by the neighbouring ranks, but no computations are performed on such rows (see the solver loop). For rank 0 (resp M_size-1) the bottom (resp. top) row is used to store the boundary condition on the bottom (resp. top) boundary. Each local matrix also stores the left and right boundary conditions (but yet those rows are only set at the beginning and never modified in the solver loop).
- In order to use MPI, it stores some methods returing the pointer to the starting position where to send/receive the top/bottom data.

`laplaceProblem`: the main class handling the problem. It implements methods to:
- initialize the problem, reading the parameters and setting the boundary conditions.
- solve the problem, performing local computations (using openMP by default), exchanging information using MPI and compute errors.
- export the solution in `vtk` format.

`MuparserFun`: the class given in the Extras of the course with some slight modifications to enable ...

## Notes for compiling the code
In the Makefile, you should change `$(PACS_ROOT)` to where your `pacs-example` folder resides.

- To compile the solver, just go in the `test` directory and run `make`.
- To build the executable, run `make test`. By default this will use a Hybrid parallelization approach.
- If you want to compare the results between a Hybrid approach and a MPI-only parallelization, run `make no_omp`. This will rebuild the executable without the `-fopenmp` flag. Notice that a comparison has been already done in `performance_hybrid.txt` and `performance_MPI.txt`.
- Run `make doc` to generate documentation.

## Optional runtime flags
To run the executable type 
```bash
mpiexec --bind-to none -n <nproc> ./main
```
with possible options:
  - `-v`: be verbose (this will drastically slow down the performances, so use it just for debug purposes).
  - `-f <filename>`: to give a differnet filename (`data.pot` by default).


To run the scalability test type
```bash
./run_scalability_test.sh
```
with possible options:
  - `-v`: be verbose (this will drastically slow down the performances, so use it just for debug purposes).
  - `-p`: visualize the results using mayavi2.

