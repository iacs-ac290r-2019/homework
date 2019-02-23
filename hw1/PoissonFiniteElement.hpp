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
    //  Constructor & Destructor
    /** Constructor: load from a YAML file
     *  @param[in] fname the name of the YAML configuration file, e.g. my_problem.yml
     */
    PoissonFiniteElement(string fname = "");

    /// Destructor; needs to delete all manually created arrays of doubles used with LAPACK
    ~PoissonFiniteElement();

    // *****************************************************************************************************************
    // Data elements
    private:
    /// The name of the configuration file used to set up the problem
    string fname;

    /// Vector of sampled values of f(x) at the node points
    vector<double> f;

    /// The boundary condition u(1) = g
    double g;

    /// The boundary condition u_x(0) = h
    double h;

    /// The number of elements, e.g. 1024
    int n;

    /// The degrees of freedom for each node, e.g. 2 for piecewise linear elements
    int d;

    /// Order of Gaussian Quadrature to use, e.g. 1 to use the midpoint
    int q;

    /** The grid of node locations; will be uniformly spaced, store a vector 
     * for extensibility to non-uniform mesh size in the future. */
    vector<double> x;

    /// The global stiffness matrix, K
    double *K;

    /// The global force vector, F
    double *F;

    // *****************************************************************************************************************
    // Accessors
    public:

    inline int num_elements() const
    {return n;}

    /** Access entry (i, j) of the global stiffness matrix K
     *  @param[in] i the row
     *  @param[in] j the column
     *  @return K_ij the stiffness K[i, j]
     */
    inline double K_ij(int i, int j) const
    {return K[ij2k(i, j)];}

    /** Access entry (i) of the global force matrix F
     *  @param[in] i the row
     *  @return F_i the force F[i]
     */
    inline double F_i(int i) const
    {return F[i];}

    // *****************************************************************************************************************
    // Calculations of indices and element size
    public:

    /** Convert a pair (i, j) of row and column indices on an mxn matrix to a storage index k 
     *  when the matrix is stored in column-major order for Fortran compatibility.
     * 
     *  @param[in] i the row
     *  @param[in] j the column
     *  @param[in] n the number of columns in the matrix
     *  @return k the storage location
     */
    inline int ij2k_calc(int i, int j, int n) const
    {return n*j + i;}

    /** Convert a pair (i, j) of row and column indices on the stiffnes matrix to a storage index k
     * 
     *  @param[in] i the row on the global stiffness matrix K
     *  @param[in] j the column on the global stiffness matrix K
     *  @return k the storage location on the global stiffness matrix K
     */
    inline int ij2k(int i, int j) const
    {return ij2k_calc(i, j, n);}

    /** Convert a pair (i, j) of row and column indices on K_element to a storage index k
     * 
     *  @param[in] i the row on the local stiffness matrix K_element
     *  @param[in] j the column on the local stiffness matrix K_element
     *  @return k the storage location on the local stiffness matrix K_element
     */
    inline int ij2k_elt(int i, int j) const
    {return ij2k_calc(i, j, d);}

    /** Get the size h of element e
     * 
     *  @param[in] e the index of the element
     *  @return h_e the size of element e; will be a constant 1.0/n for uniform mesh size
     */
    inline double h_e(int e)  const
    {return x[e+1] - x[e];}

    // *****************************************************************************************************************
    // Calculation of stiffness matrix K
    public:

    /** Create the element stiffness matrix, K_element, of size kxk, in the specific case d=2.
     * 
     *  @param[in] e the element number, e.g. 10
     *  @param[out] Ke, an array of size k^2
     * 
     *  This matrix will be stored in column major order for LAPACK compatibility.
     *  See Hughes p. 45 for details of the simple calculation.
     */
    void K_element_2(int e, double *Ke);

    /** Assemble the global stiffness matrix K as a dense nxn matrix in the general case
     *  by dispatching a call to the appropriate assembly function.*/ 
    void assemble_K();

    /** Assemble the global stiffness matrix K as a dense nxn matrix in the case d=2.*/ 
    void assemble_K_2();

    // *****************************************************************************************************************
    // Calculation of force vector F
    public:

    /** Create the element stiffness vector, F_element, of size k, in the specific case d=2
     * 
     *  @param[in] e the element number, e.g. 10
     *  @param[out] Kf, an array of size k
     */
    void F_element_2(int e, double *Fe);
    
    /** Assemble the global stiffness force vector F in the general case
     *  by dispatching a call to the appropriate assembly function.*/ 
    void assemble_F();

    /** Assemble the global force vector F in the case d=2.*/ 
    void assemble_F_2();

    // *****************************************************************************************************************
    // Output methods
    public:

    /// Write a summary of problem setup to the console
    void print_problem() const;

    /// Write a summary of stiffness matrix K
    void print_K() const;

    /// Write a summary of force vector F
    void print_F() const;

};
