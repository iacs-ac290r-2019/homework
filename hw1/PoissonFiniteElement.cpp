/**
 *  @file PoissonFiniteElement.cpp
 *  @brief Implementation of PoissonFiniteElement class.
 * 
 *  @author Michael S. Emanuel
 *  @date 2019-02-20
 */

// *********************************************************************************************************************
#include "PoissonFiniteElement.hpp"

// *********************************************************************************************************************
// Constructor - build from an input file
PoissonFiniteElement::PoissonFiniteElement(string fname) :
    fname(fname), 
    f(vector<double>()), 
    g(0.0), h(0.0), n(0), d(0), q(0),
    x(vector<double>()),
    K(nullptr)
{
    // https://stackoverflow.com/questions/45346605/example-for-yaml-cpp-0-5-3-in-linux
    try
    {  
        // Load the YAML configuration file
        YAML::Node config = YAML::LoadFile(fname);

        // Set the members to the configured values
        f = config["f"].as<vector<double>>();
        g = config["g"].as<double>();
        h = config["h"].as<double>();
        n = config["n"].as<int>();
        d = config["d"].as<int>();
        q = config["q"].as<int>();

        // Check that length of f is equal to n+1
        if (f.size() != n+1)
        {
            string msg = (format("Error: length of f %1% is not equal to n+1 = %2%.\n") % f.size() % (n+1)).str();
            throw std::runtime_error(msg);
        }
    }
    // If there was a problem above, e.g. bad file or bad input, will be left an empty PoissonFiniteElement   
    catch (const std::runtime_error &e) 
    {
        cout << e.what();
    }
    catch (...)
    {
        cout << "Unknown exception encountered in PoissonFiniteElement constructor; returning empty problem.";
    }
    // Populate the vector x of node points using a uniform mesh size
    x.resize(n+1);
    for (int i=0; i<=n; ++i)
    {
        x[i] = double(i) / n;
    }
}


PoissonFiniteElement::~PoissonFiniteElement()
{
    // Don't need to check if K was assigned memory b/c it is safe (and legal) to delete nullptr in C++17
    delete [] K;
}

// *********************************************************************************************************************
// Compute K_element; this is a 2x2 matrix for piecewise linear basis function.
// See Hughes p. 45.
void PoissonFiniteElement::K_element(int e, double *Ke)
{
    // Calculation will depend on the degrees of freedom, d
    // For now the only supported configuration is piecewise linear with d=2
    switch (d)
    {
        case 2:
            // Dispatch call to K_element_2
            K_element_2(e, Ke);
            break;
        default:
            string msg = (format("Error: degrees of frredom d=%i not supported.") % d).str();
            throw std::logic_error(msg);
    }
}


// *********************************************************************************************************************
void PoissonFiniteElement::K_element_2(int e, double *Ke)
{
    // Compute the size of this element
    double h = h_e(e);
    // Cache the inverse of h
    double h_inv = 1.0 / h;
    // The result is a symmetric 2x2 matrix; see Hughes p. 45
    // There are 1s on the main diagonal and -1s on the off diagonals
    // Everything is premultiplied by 1 / h

    // The main diagonal of a 2x2 matrix is stored on elements 0 and 3 (first and last)
    Ke[0] = h_inv;
    Ke[3] = h_inv;
    // The off diagonal of a 2x2 matrix is on the other 2 slots
    Ke[1] = -h_inv;
    Ke[2] = -h_inv;
}

// *********************************************************************************************************************
void PoissonFiniteElement::assemble_K()
{
    // Calculation will depend on the degrees of freedom, d
    // For now the only supported configuration is piecewise linear with d=2
    switch (d)
    {
        case 2:
            // Dispatch call to assemble_K_2
            assemble_K_2();
            break;
        default:
            string msg = (format("Error: degrees of frredom d=%i not supported.") % d).str();
            throw std::logic_error(msg);
    }
}

// *********************************************************************************************************************
void PoissonFiniteElement::assemble_K_2()
{
    // TODO: In the future may wish to parallelize this
    // Local storage for K_element(e)
    double Ke[4];

    // Initialize an array for K of size nxn with all zeroes
    K = new double[n*n] {0.0};

    // Iterate over the first n-1 elements, which have the same treatment
    for (int e=0; e<n-1; ++e)
    {
        // Compute K_element for this element
        K_element_2(e, Ke);
        // Increment relevant entries in global K matrix
        // Formula is on bottom of p. 42 in Hughes
        // For all Hughes formulas, shift from 1-based to 0-based indexing

        // The two diagonal entries
        // K[e,e] += K_element[1,1]
        K[ij2k(e, e)] += Ke[ij2k_elt(0, 0)];
        // K[e+1, e+1] += K_element[2, 2]
        K[ij2k(e+1, e+1)] += Ke[ij2k_elt(1, 1)];

        // The two off diagonal entries are treated symmetrically
        // K[e,e+1] ++ K_element[1, 2]
        // K[e+1,e] ++ K_element[1, 2]
        K[ij2k(e+1, e)] += Ke[ij2k_elt(1, 0)];
        K[ij2k(e, e+1)] += Ke[ij2k_elt(0, 1)];
    }

    // Handle the last element, e=n-1; only the local node (e,e) contributes
    // because node n is locked by the constraint
    int e=n-1;
    K_element_2(e, Ke);
    K[ij2k(e, e)] += Ke[ij2k_elt(0, 0)];
}

// *********************************************************************************************************************
// Print summary of problem setup
void PoissonFiniteElement::print_problem() const
{
    // Summarize this problem instance
    cout << format("Loaded 1D Poisson Problem from configuration file %1%.\n") % fname;
    cout << format("u(1)       g=%5.3f\n") % g;
    cout << format("u'(0)      h=%5.3f\n") % h;
    cout << format("mesh size  n=%i\n") % n;
    cout << format("DOF        d=%i\n") % d;
    cout << format("Quadrature q=%i\n") % q;

    // Print the f vector
    cout << format("\nf vector of %i elements:\n") % f.size();
    for (double f_e : f)
    {
        cout << format("%4.2f,  ") % f_e;
    }
    cout << "\n";

    // Print the x vector
    cout << format("\nx vector of %i nodes:\n") % x.size();
    for (double x_i : x)
    {
        cout << format("%4.2f,  ") % x_i;
    }
    cout << "\n";

}

// *********************************************************************************************************************
void PoissonFiniteElement::print_K() const
{
    cout << "\n";
    // Iterate over rows
    for (int i=0; i<n; ++i)
    {
        // Print each entry in row i
        for (int j=0; j<n; ++j)
        {
            cout << format("%5.3f ") % K_ij(i , j);
        }
        cout << "\n";
    }
}