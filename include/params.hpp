#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_
#include <string>
#include <functional>
#include <vector>
#include <array>
#include <numbers>
#include <unordered_map>
#include <string>
#include <cmath>

#include "GetPot"

namespace laplaceSolver
{
    enum Side {
        BOTTOM = 1,
        TOP = 2,
        LEFT = 3,
        RIGHT = 4
    };

    // just a forward declaration
    class MuparserFun;

    using size_type = std::size_t;
    using scalar_type = double;
    using vector_type = std::vector<scalar_type>;
    using idx_type = std::pair<size_type, size_type>;
    using node_type = std::array<scalar_type,2>;
    // using function_type = std::function<scalar_type(node_type)>;
    using BCstring_type = std::unordered_map<Side, std::string>;
    using BC_type = std::unordered_map<Side, MuparserFun>;
    constexpr auto pi = std::numbers::pi_v<scalar_type>;



    struct Numerics{

        scalar_type L = 1; ///< length of the square domain
        int n = 4; ///< number of points in each direction
        unsigned maxit = 500; ///< maximum number of iterations
        scalar_type tol = 1.e-4; ///< solver tolerance
        bool verbose = false; ///< verbosity flag

    };

    /**
     * @brief Struct that holds all the parameters for the problem
     * Parameters can be either read from a file using GetPot otherwise their are defaulted
     */
    struct Params
    {
        
        Numerics np;

        std::string f = "8*pi*pi*sin(2*pi*x1)*sin(2*pi*x2)";
        std::string u_bottom;
        std::string u_left;
        std::string u_right;
        std::string u_top;
        std::string sol_exact = "sin(2*pi*x1)*sin(2*pi*x2)";

    };

   

    /**
     * @brief Overload of the output stream operator
     */
    std::ostream& operator<<(std::ostream& os, const Params&);

    /**
     * @brief Reads problem parameters from GetPot file
     * @param filename The getopot file with the new values
     * @param verbose Prints some information on the parameters
     * If the datafile is not found, default values for the parameters are used
     */
    Params readParameters(const std::string &filename, bool verbose);

} // namespace laplaceSolver


#endif // _PARAMS_HPP_