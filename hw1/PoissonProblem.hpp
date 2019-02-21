/**
 *  @file PoissonProblem.hpp
 *  @brief Class defining one instance of the 1D Poisson problem, i.e. boundary values, f(x), and
 *  the number of grid points.
 * 
 *  @author Michael S. Emanuel
 *  @date 2019-02-20
 */

// *********************************************************************************************************************
#pragma once

// Silence selected GCC warnings
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wreorder"
#endif

// *********************************************************************************************************************
// Libraries
#include <vector>
    using std::vector;

#include <string>
    using std::string;

// YAML
#include<yaml-cpp/yaml.h>

// *********************************************************************************************************************
/**
 * @class PoissonProblem
 *  @brief Class defining one instance of the 1D Poisson problem, i.e. boundary values, f(x), and
 *  the number of grid points.
 */
class PoissonProblem
{
    public:
        /** Constructor: load from a YAML file
         * @param[in] fname the name of the YAML configuration file, e.g. my_problem.yml
         */
        PoissonProblem(string fname);

        /// Vector of sampled values of f(x) at the node points
        vector<double> f;

        /// The boundary condition u(1) = g
        double g;

        /// The boundary condition u_x(0) = h
        double h;

        /// The number of elements, e.g. 1024
        int n;

        /// The degrees of freedom for each node, e.g. 2 for piecewise linear elements
        int k;

        /// Order of Gaussian Quadrature to use, e.g. 1 to use the midpoint
        int q;

    private:

};
