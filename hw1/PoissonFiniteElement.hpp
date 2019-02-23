/**
 *  @file PoissonFiniteElement.hpp
 *  @brief Solve the 1D Poisson Problem using the Finite Element Method.
 * 
 *  Poisson Problem is 
 *  (1) u_xx + f(x) = 0
 *  (2) -u_x(0) = h
 *  (3) u(1) = g
 *  Follows the treatment of "The Finite Element Method", Thomas Hughes, chapter 1.
 * 
 *  @author Michael S. Emanuel
 *  @date 2019-02-20
 */

// *********************************************************************************************************************
#pragma once

// Silence selected GCC warnings
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wcatch-value="
#endif

// *********************************************************************************************************************
// Libraries
#include <vector>
    using std::vector;

#include <string>
    using std::string;

#include <iostream>
    using std::cout;

#include <boost/format.hpp>
    using boost::format;

// YAML
#include <yaml-cpp/yaml.h>

// *********************************************************************************************************************
/**
 * @class PoissonFiniteElement
 * @brief Class to define and solve one instance of the 1D Poisson problem.
 */
class PoissonFiniteElement
{
    public:
    // *****************************************************************************************************************
    // Constructor & destructor
    /** Default constructor: create an empty problem instance     */
    // PoissonFiniteElement();

    /** Constructor: load from a YAML file
     * @param[in] fname the name of the YAML configuration file, e.g. my_problem.yml
     */
    PoissonFiniteElement(string fname = "");

    // *****************************************************************************************************************
    // Calculations
    /** Create the element stiffness matrix, K_element, of size kxk
     * 
     * @param[in] i the element row number, e.g. 1 or 2
     * @param[in] j the element column number, e.g. 1 or 2
     * @param[out] Ke, an array of size k^2
     * 
     * This matrix will be stored in column major order for LAPACK compatibility.
     */
    void K_element(int i, int j, double *Ke);

    // *****************************************************************************************************************
    // Output methods
    /// Write a summary to the console
    void summary() const;

    // *****************************************************************************************************************
    // Data elements
    private:
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

    /// The name of the configuration file
    string fname;

};
