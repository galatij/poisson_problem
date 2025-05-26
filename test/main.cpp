// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wsuggest-override"
// #pragma GCC diagnostic ignored "-Wcast-function-type"
// #pragma GCC diagnostic pop

#include <iostream>
#include <string>
#include <mpi.h>
#include "chrono.hpp"
#include "laplace.hpp"
#include <iomanip>

int main(int argc, char * argv[])
{
    using namespace laplaceSolver;

    Timings::Chrono clock;

    MPI_Init(&argc, &argv);

    laplaceSolver::laplaceProblem solver(MPI_COMM_WORLD);

    // Read params, broadcast and setup the system
    solver.init(argc, argv);

    // Solve using Jacobi iteration
    clock.start();
    solver.solve();
    clock.stop();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
        std::cout << "Solver time: " << clock.wallTime()*1.E-6
            << std::fixed << std::setprecision(9)
            << " seconds" << std::endl;

    // Write results for post-process in Paraview
    solver.exportVtk();

    MPI_Finalize();

    return 0;
}



