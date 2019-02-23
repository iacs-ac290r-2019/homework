/**
 *  @file poisson.cpp
 *  @project poisson
 *  @brief Solve the 1D Poisson problem using the Finite Element method for AC 290 Homework 1, Problem 4.
 * 
 *  @author Michael S. Emanuel
 *  @date 2019-02-20
 */

// *********************************************************************************************************************
// Libraries
#include <iostream>
    using std::cout;

#include <string>
    using std::string;

#include <boost/format.hpp>
    using boost::format;

// Local dependencies
#include "PoissonFiniteElement.hpp"

// *********************************************************************************************************************
/** Run a test case for the 1D Poisson problem
 *  @param[in] fname the name of the configuration file
 *  @param[in] U_func function with the known solution
 *  @return isOK true for pass, false for fail
 */ 
bool test(string fname, double (*U_func) (double), bool verbose = false)
{
    // Status update
    cout << format("\nLoading configuration file %1%:\n") % fname;    

    // Set up the problem instance from the configuration file
    // Need a try / catch block because bad input throws a runtime error
    PoissonFiniteElement pfe{PoissonFiniteElement(fname)};

    // Print summary of problem setup to screen
    if (verbose) {pfe.print_problem();}

    // Assemble the global stiffness matrix K
    pfe.assemble_K();

    // Assemble the global force vector F
    pfe.assemble_F();

    // Display the stiffness matrix K
    if (verbose)
    {
        cout << format("Assembled %1% x %1% stiffness matrix K:") % pfe.num_elements();
        pfe.print_K();
    }

    // Display the force vector F
    if (verbose)
    {
        cout << format("\nAssembled %1% x 1 force vector F:") % pfe.num_elements();
        pfe.print_F();
    }

    // Solve the system
    pfe.solve();

    // Display the solution u(x)
    cout << format("\nAssembled %1% x 1 solution vector U:") % pfe.num_elements();
    pfe.print_U();

    // Compare this to known solution
    int n = pfe.num_elements();
    double max_err = 0.0;
    for (int i=0; i <= n; ++ i)
    {
        // The x value at this point
        double x = double(i) / n;
        // The correct answer for U(x) at this point
        double ans = (*U_func)(x);
        // The calculatd answer at this point
        double calc = pfe.U_i(i);
        // The error at this point
        double err = std::abs(calc - ans);
        // The max error seen
        max_err = err > max_err ? err : max_err;
        // Status
        // cout << format("i=%i, x=%5.3f, calc=%5.3f, ans=%5.3f, err=%.4e\n") % i % x % calc % ans % err;
    }

    // Tolerance for passing
    double tol = 1.0E-4;
    bool isOK = max_err < tol;

    // Report max error and pass / fail
    string msg = isOK ? "PASS" : "FAIL";
    cout << format("\nMaximum error is %.4e\n") % max_err;
    cout << format("**** %s ****\n") % msg;

    // Normal program exit
    return isOK;
}


// *********************************************************************************************************************
// Test case 1 with f=0, g=0, h=1; solution is u(x) = 1-x
double U_func_1(double x)
{
    return 1.0 - x;
}

// *********************************************************************************************************************
// Test case 2 with f=0, g=1, h=1; solution is u(x) = 2-x
double U_func_2(double x)
{
    return (2.0 - x);
}

// *********************************************************************************************************************
// Test case 3 with f=1, g=1, h=1; solution is u(x) = 2-x + 1/2(1-x^2)
double U_func_3(double x)
{
    return 2.0 - x + 0.5 * (1 - x * x);
}

// *********************************************************************************************************************
void run_tests()
{
    // Test case 1 with f=0, g=0, h=1; solution is u(x) = 1-x
    test("example_1.yml", &U_func_1);

    // Test case 2 with f=0, g=1, h=1; solution is u(x) = 2-x
    test("example_2.yml", &U_func_2);

    // Test case 3 with f=1, g=1, h=1; solution is u(x) = 2-x + 1/2(1-x^2)
    test("example_3.yml", &U_func_3);
}

// *********************************************************************************************************************
int main()
{
    // Run test suite
    run_tests();

    // Name of the configuration file
    // TODO: change this to command line argument
    // string fname("example_1.yml");
    // Status update
    // cout << format("Loading configuration file %1%:\n") % fname;    

    // Normal program exit
    return 0;
}
