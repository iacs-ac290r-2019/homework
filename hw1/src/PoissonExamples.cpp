/**
 *  @file PoissonExamples.cpp
 *  @project poisson
*  @brief Force functions and solutions to simple examples of the Poisson equation
  * 
 *  @author Michael S. Emanuel
 *  @date 2019-02-23
 */

// *********************************************************************************************************************
#include "PoissonExamples.hpp"

// *********************************************************************************************************************
// Test case 1 with f=0, g=0, h=1; solution is u(x) = 1-x
double F_func_1(double x)
{
    return 0.0;
}

/// Solution function for test case 1
double U_func_1(double x)
{
    return 1.0 - x;
}

// *********************************************************************************************************************
// Test case 2 with f=0, g=1, h=1; solution is u(x) = 2-x
double F_func_2(double x)
{
    return 0.0;
}

double U_func_2(double x)
{
    return (2.0 - x);
}

// *********************************************************************************************************************
// Test case 3 with f=1, g=1, h=1; solution is u(x) = 2-x + 1/2(1-x^2)
double F_func_3(double x)
{
    return 1.0;
}

double U_func_3(double x)
{
    return 2.0 - x + 0.5 * (1 - x * x);
}

// *********************************************************************************************************************
// Test case 4; solution is u(x) -->
// u''(x) = -sin(x) --> f(x) = sin(x)
//  h=-u'(0)=-1, g=u(1)=sin(1)=0.841470984808
double F_func_4(double x)
{
    return sin(x);
}

double U_func_4(double x)
{
    return sin(x);
}

// *********************************************************************************************************************
// Factory function to return a function pointer given name of a force function

// Create a map keyed by function name, with value a pointer to the function
map<string, FuncType> funcTable {};

// Force functions
void makeFuncTable()
{
    funcTable["F_func_1"] = &F_func_1;
    funcTable["F_func_2"] = &F_func_2;
    funcTable["F_func_3"] = &F_func_3;
    funcTable["F_func_4"] = &F_func_4;

    // Solution functions
    funcTable["U_func_1"] = &U_func_1;
    funcTable["U_func_2"] = &U_func_2;
    funcTable["U_func_3"] = &U_func_3;
    funcTable["U_func_4"] = &U_func_4;
}

// Wrap this into a function
FuncType FuncByName(string funcName)
{
    return funcTable[funcName];
}
