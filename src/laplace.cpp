#include "laplace.hpp"
#include "muparser_fun.hpp"
#include <numeric>
#include "mpi_utils.hpp"

bool DEBUGgath = true;

namespace laplaceSolver {

    void print_idx_vector(const std::vector<idx_type>& idxList);

    void
    laplaceProblem::init(int argc, char* argv[])
    {
        bool verbose = false;

        // read the parameters
        if (M_rank == 0) {  // rank 0 does the job
            GetPot cl(argc, argv);
            // check if we want verbosity
            verbose = cl.search("-v");// true if -v is found
            // get file with parameter values
            std::string filename = cl.follow("data.pot",2,"--file", "-f");// "data.pot" is the default, otherwise what follows -f
            M_params = readParameters(filename, verbose);
            
            if (verbose)
                std::cout << "Reading parameters...\n" << M_params << std::endl;

        }

        // broadcast the numerical parameters to all ranks
        MPI_Bcast(&M_params.np, sizeof(Numerics), MPI_BYTE, 0, M_comm);

        // concatenated string to be sent
        std::string allFunctions;

        // size of each function string
        std::vector<int> allFuncSize (6);

        if (M_rank == 0){
            // concatenate function strings
            allFunctions = M_params.f + M_params.u_bottom + M_params.u_right
                + M_params.u_top + M_params.u_left + M_params.sol_exact;

            // save the size of each function string
            allFuncSize[0] = M_params.f.size();
            allFuncSize[1] = M_params.u_bottom.size();
            allFuncSize[2] = M_params.u_right.size();
            allFuncSize[3] = M_params.u_top.size();
            allFuncSize[4] = M_params.u_left.size();
            allFuncSize[5] = M_params.sol_exact.size();

        }

        // Broadcast the sizes of the function strings
        MPI_Bcast(allFuncSize.data(),6,MPI_INT,0,M_comm);
        int totalSize = std::accumulate(allFuncSize.begin(), allFuncSize.end(),0);
        allFunctions.resize(totalSize);
        
        // Broadcast the concatenated string
        MPI_Bcast(allFunctions.data(),totalSize,MPI_CHAR,0,M_comm);

        if (M_rank != 0){
            // save the strings in M_params
            int offset = 0;
            M_params.f = allFunctions.substr(0, allFuncSize[0]);
            offset += allFuncSize[0];
            M_params.u_bottom = allFunctions.substr(offset, allFuncSize[1]);
            offset += allFuncSize[1];
            M_params.u_right = allFunctions.substr(offset, allFuncSize[2]);
            offset += allFuncSize[2];
            M_params.u_top = allFunctions.substr(offset, allFuncSize[3]);
            offset += allFuncSize[3];
            M_params.u_left = allFunctions.substr(offset, allFuncSize[4]);
            offset += allFuncSize[4];
            M_params.sol_exact = allFunctions.substr(offset, allFuncSize[5]);

        }

        // parse f and BC using muParser
        M_f.SetExpression(M_params.f);
        M_sol_exact.SetExpression(M_params.sol_exact);
        if (M_rank == 0)
            M_BC[BOTTOM].SetExpression(M_params.u_bottom);
        else if (M_rank == M_size-1)
            M_BC[TOP].SetExpression(M_params.u_top);
        
        M_BC[LEFT].SetExpression(M_params.u_left);
        M_BC[RIGHT].SetExpression(M_params.u_right);

        // setup the matrices
        std::vector<int> counts_send;
        std::vector<int> displacements;
        if (M_rank ==0) {// root has to prepare some data
            counts_send.resize(M_size);
            displacements.resize(M_size, 0);// init. by zero
            int chunk = (M_params.np.n - 2) / M_size; // - 2 to account for first and last row (BC)
            int rest = (M_params.np.n - 2) % M_size;
            for (int i{}; i < M_size;++i){ // count the nrows in each rank
                counts_send[i] = i < rest ? chunk + 1 : chunk;
                if(i > 0)
                    displacements[i] = displacements[i-1] + counts_send[i-1];
            }
        }

        int local_nrows; // each process gets local size with normal scatter
        MPI_Scatter(counts_send.data(),1,MPI_INT, &local_nrows,1,MPI_INT,0,M_comm);

        int start_row;
        MPI_Scatter(displacements.data(), 1, MPI_INT, &start_row, 1, MPI_INT, 0, M_comm);

        // resize the matrix and precompute boundary (or ghost) indices
        M_localU.initResize(local_nrows + 2, M_params.np.n, start_row); // +2 to account for ghost rows

        // assign DirBC
        if (verbose && M_rank == 0) std::cout << "Setting boundary conditions... ";
        
        setBC();
        
        if (verbose && M_rank == 0) std::cout << "done." << std::endl;
        
    }

    void
    laplaceProblem::setBC()
    {
        // set top and bottom boundary condition
        if (M_rank == 0) {
            std::vector<idx_type> bottom_idx = M_localU.getIdxList<Side::BOTTOM>(); // take bottom idx
            for (const auto& [i,j]: bottom_idx) {
                node_type x = M_localU.getCoords(i, j, M_params.np.L); // take points for each idx
                M_localU(i,j) = M_BC[BOTTOM](x); // impose Dirichlet BC
            }
        }

        else if (M_rank == M_size- 1) {
            std::vector<idx_type> top_idx = M_localU.getIdxList<Side::TOP>(); // take top idx
            for (const auto& [i,j]: top_idx) {
                node_type x = M_localU.getCoords(i, j, M_params.np.L); // take points for each idx
                M_localU(i,j) = M_BC[TOP](x); // impose Dirichlet BC
            }
        }

        // set left and right boundary conditions
        std::vector<idx_type> left_idx = M_localU.getIdxList<Side::LEFT>(); // take left idx
        for (auto it = left_idx.begin() + 1; it != left_idx.end() - 1; ++it) {
            const auto& [i, j] = *it;
            node_type x = M_localU.getCoords(i, j, M_params.np.L);
            M_localU(i,j) = M_BC[LEFT](x);
        }

        std::vector<idx_type> right_idx = M_localU.getIdxList<Side::RIGHT>(); // take right idx
        for (auto it = right_idx.begin() + 1; it != right_idx.end() - 1; ++it) {
            const auto& [i, j] = *it;
            node_type x = M_localU.getCoords(i, j, M_params.np.L);
            M_localU(i,j) = M_BC[RIGHT](x);
        }

    }


    void
    laplaceProblem::solve()
    {
        // set number of threads to be used for each rank
        #ifdef _OPENMP
            omp_set_num_threads(omp_get_num_procs()/M_size);
        #pragma omp parallel
            #pragma omp master
            if (M_rank == 0)
                std::cout << "Number of threads for each rank: " << omp_get_num_threads() << std::endl;
        #endif
        
        if (M_rank == 0)
            std::cout << "Number of points: " << M_localU.ncols()
                << "\nSolving with Jacobi iteration...\n";

        // 1. Initialize the solution
        // copy construct the built matrix
        SolMatrix U0 {M_localU};

        size_type nrows = M_localU.nrows();
        size_type ncols = M_localU.ncols();
        
        scalar_type h = M_params.np.L/(ncols-1);
        int converged = 0;
        size_type k{};
        scalar_type totalRes{};

        // precompute forcing term
        std::vector<scalar_type> fun;
        for (size_type i = 0; i < nrows; ++i)
            for (size_type j = 0; j < ncols ; ++j)
                fun.emplace_back(M_f(M_localU.getCoords(i,j, M_params.np.L)));
       
        for ( ; k < M_params.np.maxit && !converged; ++k) {
            if (M_params.np.verbose && M_rank == 0)
                std::cout << "-- it = " << k;

            // 2. Perform local computations
            #pragma omp parallel for
            for (size_type i = 1; i < nrows - 1; ++i)
                for (size_type j = 1; j < ncols - 1; ++j)
                    M_localU(i,j) = 0.25*( U0(i-1,j) + U0(i+1,j) + U0(i,j+1) + U0(i,j-1) + 
                                         h*h*fun[i*ncols+j]);

            // 3. Send the computed values at the border to the adjacent ranks and receive from them
            exchangeGhostRows();

            // 4. Compute the local error and check local convergence criterion
                /* E.g.: each rank computes on its internal + boundary nodes*/
            scalar_type err{};
            size_type jend = ( (M_rank == M_size-1) ? nrows : nrows-1);

            #pragma omp parallel for reduction(+:err) 
            for (size_type j = (M_rank == 0 ? 0 : 1); j < jend; ++j)
                for (size_type i{}; i < ncols; ++i)
                    err += (M_localU(j,i) - U0(j,i))*(M_localU(j,i) - U0(j,i));
            
            err = std::sqrt(h*h*err); // Using h^2 (differently from pdf)

            // check convergence criterion
            converged = (err < M_params.np.tol*h) ? 1 : 0;
            
            if (M_params.np.verbose){
                /** @remark Printing the residual at each iteration requires
                 * more communication, slowing down the solver
                 */
                scalar_type err2 = err*err;
                MPI_Reduce(&err2, &totalRes, 1, mpi_typeof(err2), MPI_SUM, 0, M_comm);
                if (M_rank == 0)
                    std::cout << "\t ||res|| = " << std::sqrt(totalRes) << std::endl;
            }

            MPI_Allreduce(MPI_IN_PLACE, &converged, 1, MPI_INT, MPI_PROD, M_comm);

            // 6. Update solution
            U0 = M_localU;
        }

        // 7. Print information
        if (converged){
            M_hasConverged = true;
            scalar_type totalError{};
            SolMatrix U{this->buildGlobalMatrix(),ncols,ncols};

            if (M_rank == 0)
                totalError = errorFromExact(U, M_sol_exact);

            if (M_rank == 0){
                std::cout << "Converged in " << k << " iterations! STONKS" << std::endl;
                if (M_params.np.verbose)
                    std::cout << "Final residual: " << std::sqrt(totalRes) 
                        << ", Error: " << totalError
                        << std::endl;
            }
        }
        else if (M_rank == 0)
            std::cout << "Solver DIDN'T converge" << std::endl;

    }


    scalar_type laplaceProblem::errorFromExact(SolMatrix U_sol, MuparserFun sol_exact){
        /** \todo: Try to parallelize this avoiding datarace in muParser */
        //  #ifdef _OPENMP
        //  int n_threads=omp_get_num_procs()/M_size; // così n=2 non è ottimizzato
            //  if(M_size==2)
            //     n_threads=2;
        //     #endif 
        scalar_type err{};
        //# pragma omp parallel for num_threads(n_threads) reduction(+:err)
        for (size_type i=0; i < U_sol.nrows();++i){
            for (size_type j=0; j < U_sol.ncols();++j){
                scalar_type val = (U_sol(i,j) - sol_exact(U_sol.getCoords(i,j,M_params.np.L)));// Problemi con muparser, vedere se non viene solo letto
                err += val*val;
            }
        }
    
        scalar_type h = M_params.np.L/(U_sol.ncols()-1);
        return err = std::sqrt(h*h*err); // Using h^2 (differently from slides)

    }


    void laplaceProblem::exchangeGhostRows(){
        const auto &nsub = M_params.np.n;
        if (M_rank != 0) { // send/receive bottom row
            int adjacentRank = M_rank - 1;
            MPI_Sendrecv(M_localU.sendBottomData(), nsub-2, mpi_typeof(M_localU(0,0)), adjacentRank, 1,
                         M_localU.recvBottomData(), nsub-2, mpi_typeof(M_localU(0,0)), adjacentRank, 2,
                         M_comm, MPI_STATUS_IGNORE);
        }
        if (M_rank != M_size-1) { // send/receive top row
            int adjacentRank = M_rank + 1;
            MPI_Sendrecv(M_localU.sendTopData(), nsub-2, mpi_typeof(M_localU(0,0)), adjacentRank, 2,
                         M_localU.recvTopData(), nsub-2, mpi_typeof(M_localU(0,0)), adjacentRank, 1,
                         M_comm, MPI_STATUS_IGNORE);
        }
    }

    /** @remark only for DEBUG: */
    void print_idx_vector(const std::vector<idx_type>& idxList){
        for (const auto& [i,j]: idxList)
            std::cout << "(" << i << ","<< j <<")" << std::endl;
    }



    std::vector<scalar_type>
    laplaceProblem::buildGlobalMatrix() const
    {
        // gather the results on rank 0 for output
        size_type dim = M_params.np.n;
        std::vector<scalar_type> globalU;

        std::vector<int> counts_recv, displacements;

        if (M_rank == 0) { // root has to prepare some data
            counts_recv.resize(M_size);
            displacements.resize(M_size);

            globalU.resize(dim*dim);
            int chunk = (dim - 2) / M_size; // - 2 to account for first and last row (BC)
            int rest = (dim - 2) % M_size;
            for (int i{}; i < M_size; ++i){ // count the number of elements to be sent in each rank
                counts_recv[i] = i < rest ? dim*(chunk + 1) : dim*chunk;
                if (i == 0 || i == M_size -1)
                    counts_recv[i] += dim; // the first/last rank sends also the bottom/top boundary
                if(i > 0)
                    displacements[i] = displacements[i-1] + counts_recv[i-1];
            }
        }

        int sendSize = M_localU.effectiveSize();
        sendSize += (M_rank == 0 || M_rank == M_size-1) ? dim : 0;

        // Collect local computation to rank 0
        MPI_Gatherv(M_localU.data(M_rank == 0), 
                sendSize,
                mpi_typeof(M_localU(0,0)),
                M_rank == 0 ? globalU.data() : nullptr,
                counts_recv.data(),
                displacements.data(),
                mpi_typeof(M_localU(0,0)),
                0,
                M_comm);

        return globalU;
    }
    
    void
    laplaceProblem::exportVtk(const std::string& filename) const
    {
        // gather the local matrix in rank 0
        size_type dim = M_params.np.n;
        auto U = buildGlobalMatrix();
        
        if (M_rank == 0){
            // open the file for writing
            std::ofstream vtkFile(filename);
            if (!vtkFile)
                throw std::runtime_error("Cannot open output file: " + filename);

            vtkFile << "# vtk DataFile Version 3.0\n";
            vtkFile << "Matrix output\n";
            vtkFile << "ASCII\n";
            vtkFile << "DATASET STRUCTURED_POINTS\n";
            vtkFile << "DIMENSIONS " << dim << " " << dim << " 1\n";
            vtkFile << "ORIGIN 0 0 0\n";
            vtkFile << "SPACING 1 1 1\n";
            vtkFile << "POINT_DATA " << dim * dim << "\n";
            vtkFile << "SCALARS " << "u" << " double\n";
            vtkFile << "LOOKUP_TABLE default\n";

            // write data in row-major order
            for (size_t i {}; i < dim; ++i)
                for (size_t j {}; j < dim; ++j)
                    vtkFile << U[i * dim + j] << "\n";

            vtkFile.close();
        }
    }


} // namespace laplaceProblem