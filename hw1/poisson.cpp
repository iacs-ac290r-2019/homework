/**
 *  @file poisson.cpp
 *  @project poisson
 *  @brief Solve the 1D Poisson problem using the Finite Element method for AC 290 Homework 1, Problem 4.
 * 
 *  @author Michael S. Emanuel
 *  @date 2019-02-20
 */

// *********************************************************************************************************************
#include <iostream>
    using std::cout;

#include <boost/format.hpp>
    using boost::format;

#include<yaml-cpp/yaml.h>

// *********************************************************************************************************************
int main()
{
    int n =42;
    cout << format("Hello.  n=%1%\n") % n;
    return 0;
}
