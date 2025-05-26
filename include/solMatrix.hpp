#ifndef _SOL_MATRIX_HPP_
#define _SOL_MATRIX_HPP_

#include "params.hpp"

namespace laplaceSolver {

    /**
     * @brief A Matrix class to store the solution in each grid node
     * It has specific methods to enable parallelization and exchange of information
     * for the solution of the laplace problem
     */
    class SolMatrix {
        size_type M_nrowsInternal;
        size_type M_nrowsWithGhost;
        size_type M_globalStartRow;
        size_type M_ncols;

        std::vector<scalar_type> M_mat;

    public:

        SolMatrix(std::vector<scalar_type> v, size_type nrows,size_type ncols)
        : M_nrowsInternal(nrows-2),
        M_nrowsWithGhost(nrows),
        M_globalStartRow(0),
        M_ncols(ncols),
        M_mat(v)
        {}

        SolMatrix() = default;

        /**
         * @brief Prepares the matrix with correct sizes
         */
        void initResize(size_type nrowsWithGhost, size_type ncols, size_type start_row){
            M_globalStartRow = start_row;
            M_nrowsInternal = nrowsWithGhost-2;
            M_nrowsWithGhost = nrowsWithGhost;
            M_ncols = ncols;
            M_mat.resize(nrowsWithGhost*ncols);
        }

        /** 
         * @brief Overload of the call operator
         * @param i The local row index
         * @param j The local col index
         * @return A(i,j) of the local matrix
         */
        scalar_type operator()(size_type i, size_type j) const {
            return M_mat[i*M_ncols + j];
        }

        scalar_type& operator()(size_type i, size_type j) {
             return M_mat[i*M_ncols + j];
        }

        size_type globalRow(size_type local_i) const {
            return M_globalStartRow + local_i; // assuming ghost at index 0
        }

        idx_type globalIndices(size_type local_i, size_type j) const {
            return {globalRow(local_i), j};
        }

        node_type getCoords(size_type i, size_type j, scalar_type L) const {
            scalar_type h = L/(M_ncols-1);
            scalar_type x = j*h;
            scalar_type y = globalRow(i)*h;
            return node_type{x,y};
        }
        
        scalar_type* sendBottomData() {
            return &M_mat[M_ncols+1];
        }

        scalar_type* recvBottomData() {
            return &M_mat[1];
        }

        scalar_type* sendTopData() {
            // return &M_mat[(M_nrowsWithGhost - 2) * M_ncols];
            return &M_mat[(M_nrowsWithGhost - 2)*M_ncols + 1];
        }

        scalar_type* recvTopData() {
            return &M_mat[(M_nrowsWithGhost - 1)*M_ncols + 1];
        }

        template<Side S>
        std::vector<idx_type> getIdxList() const {
            std::vector<idx_type> idx_vector;
            if constexpr (S == Side::BOTTOM) {
                // return std::vector<idx_type>(M_mat.cbegin(), M_mat.cbegin()+M_ncols);
                idx_vector.reserve(M_ncols);
                for (size_type j{}; j < M_ncols; ++j)
                        idx_vector.emplace_back(0, j);
            }

            else if constexpr (S == Side::TOP) {
                // return std::vector<idx_type>(M_mat.crbegin(), M_mat.crbegin()+M_ncols);
                idx_vector.reserve(M_ncols);
                for (size_type j{}; j < M_ncols; ++j)
                        idx_vector.emplace_back(M_nrowsWithGhost-1, j);
            }
            
            else if constexpr (S == Side::LEFT) {
                idx_vector.reserve(M_nrowsWithGhost);
                for (size_type i{}; i < M_nrowsWithGhost; ++i)
                    idx_vector.emplace_back(i,0);
            }

            else if constexpr (S == Side::RIGHT) {
                idx_vector.reserve(M_nrowsWithGhost);
                for (size_type i{}; i < M_nrowsWithGhost; ++i)
                    idx_vector.emplace_back(i,M_ncols-1);
            }

            return idx_vector;
            
        }

        scalar_type* data(bool isRank0 = false) {
            if (isRank0)
                return M_mat.data();
            return &M_mat[M_ncols];
        }

        const scalar_type* data(bool isRank0 = false) const {
            if (isRank0)
                return M_mat.data();
            return &M_mat[M_ncols];
        }

        size_type nrows() const { return M_nrowsWithGhost; }

        size_type ncols() const { return M_ncols; }

        size_type effectiveSize() const {
            return M_ncols*M_nrowsInternal;
        };

    };

} // namespace laplaceSolver


#endif // _SOL_MATRIX_HPP_


