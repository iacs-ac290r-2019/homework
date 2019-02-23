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

#include <boost/format.hpp>
    using boost::format;

// Local dependencies
#include "PoissonFiniteElement.hpp"

// *********************************************************************************************************************
int main()
{
    // Name of the configuration file
    // TODO: change this to command line argument
    string fname("example_1.yml");
    // Status update
    cout << format("Loading configuration file %1%:\n") % fname;    

    // Set up the problem instance from the configuration file
    // Need a try / catch block because bad input throws a runtime error
    PoissonFiniteElement pfe{PoissonFiniteElement(fname)};

    // Print summary of problem setup to screen
    pfe.print_problem();

    // Assemble the global stiffness matrix K
    pfe.assemble_K();

    // Assemble the global force vector F
    pfe.assemble_F();

    // Display the stiffness matrix K
    cout << format("\nAssembled %1% x %1% stiffness matrix K:") % pfe.num_elements();
    pfe.print_K();

    // Display the force vector F
    cout << format("\nAssembled %1% x 1 force vector F:") % pfe.num_elements();
    pfe.print_F();

    // Normal program exit
    return 0;
}
