#ifndef _LAPLACE_HPP_
#define _LAPLACE_HPP_


#include "params.hpp"
#include "solMatrix.hpp"
#include "muparser_fun.hpp"
#include "GetPot"
#include <mpi.h>
#include <omp.h>

namespace laplaceSolver {

    /**
     * @brief A class to solve the laplace problem in parallel
     */
    class laplaceProblem {
        MPI_Comm M_comm;
        int M_rank;
        int M_size;

        Params M_params; ///< The problem parameters
        BC_type M_BC; ///< A map of BC functors (MuParserFun)
        MuparserFun M_f; ///< The forcing term
        MuparserFun M_sol_exact; ///< The exact solution
        SolMatrix M_localU; ///< The matrix to store the (local) computed solution
        bool M_hasConverged = false; ///< flag

        /**
         * @brief A function to set the Boundary conditions in the initialization phase
         * Each rank sets the BC of of the nodes of which it is resposible
         */
        void setBC();

        /**
         * @brief Called after local computations
         */
        void exchangeGhostRows();

        /**
         * @brief Collect the local matrices to a global matrix in rank 0
         * If the solver converges, rank 0 gathers the local matrices
         * in order to allow exporting for visualization
         */
        std::vector<scalar_type> buildGlobalMatrix() const;

    public:

        laplaceProblem(MPI_Comm comm)
        : M_comm(comm) {
            MPI_Comm_rank(comm, &M_rank);
            MPI_Comm_size(comm, &M_size);
        };

        /**
         * @brief Initalizes the problem
         * Reads the parameters and sets boundary conditions
         */
        void init(int argc, char* argv[]);

        /**
         * @brief Solve the problem using Jacobi iteration
         */
        void solve();

        bool converged() const { return M_hasConverged; }

        /**
         * @brief Export the solution for vtk visuualization
         */
        void exportVtk(const std::string& filename = "output.vtk") const;

        /**
         * @brief Computes the error
         * Called at the end of the the solver algorithm by rank 0
         */
        scalar_type errorFromExact(SolMatrix U_sol, MuparserFun sol_exact);
    };

} // namespace laplaceSolver

#endif // _LAPLACE_HPP_
