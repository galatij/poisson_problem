#include "params.hpp"
#include <fstream>


namespace laplaceSolver {


    std::ostream&
    operator<<(std::ostream &out, const Params &p)
    {
        out << "- L: " << p.np.L
            << "\n- Number of subdivisions: "<< p.np.n
            << "\n- maxit: " << p.np.maxit
            << "\n- tol: " << p.np.tol
            << "\n- force: " << p.f
            << "\n- u_left: " << p.u_left
            << "\n- u_right: " << p.u_right
            << "\n- u_top: " << p.u_top
            << "\n- u_bottom: " << p.u_bottom;
        return out;
    }

    Params
    readParameters(const std::string &filename, bool verbose)
    {
        // Parameter default constructor fills it with the defaults values
        Params defaults;
        // checks if file exixts and is readable
        std::ifstream check(filename);
        if(!check) {
            std::cerr << "ERROR: Parameter file " << filename << " does not exist"
                        << std::endl;
            std::cerr << "Reverting to default values." << std::endl;
            check.close();
            return defaults;
        }
        else
            check.close();
        
        GetPot ifile(filename.c_str());
        
        Params values;

        // Read parameters from getpot data base
        values.np.L = ifile("L", defaults.np.L);
        values.np.n = ifile("n", defaults.np.n);
        values.np.tol = ifile("tol", defaults.np.tol);
        values.np.maxit = ifile("maxit", defaults.np.maxit);
        values.f = ifile("f", defaults.f.c_str());
        values.u_bottom = ifile("u1", defaults.u_bottom.c_str());
        values.u_right = ifile("u2", defaults.u_right.c_str());
        values.u_top = ifile("u3", defaults.u_top.c_str());
        values.u_left= ifile("u4", defaults.u_left.c_str());

        values.np.verbose = verbose;
        
        return values;
    }

} // namespace laplaceSolver

