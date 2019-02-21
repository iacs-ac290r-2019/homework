/**
 *  @file PoissonProblem.cpp
 *  @brief Implementation of PoissonProblem class.
 * 
 *  @author Michael S. Emanuel
 *  @date 2019-02-20
 */

// *********************************************************************************************************************
#include "PoissonProblem.hpp"

// *********************************************************************************************************************
// Constructor
PoissonProblem::PoissonProblem(string fname) :
    f(vector<double>()), g(0.0), h(0.0), n(0), k(0), q(0)
{
    // Load the YAML configuration file
    // TODO: replace hard coded example file name with fname
    YAML::Node config = YAML::LoadFile("example_1.yml");

    // Set the members to the configured values

}