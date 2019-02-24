/**
 *  @file PoissonExamples.hpp
 *  @project poisson
 *  @brief Force functions and solutions to simple examples of the Poisson equation
 * 
 *  @author Michael S. Emanuel
 *  @date 2019-02-23
 */

// *********************************************************************************************************************
#pragma once

// *********************************************************************************************************************
#include <cmath>

#include <string>
    using std::string;

#include <map>
    using std::map;

// *********************************************************************************************************************
// Test case 1 with f=0, g=0, h=1; solution is u(x) = 1-x
/// Force function for test case 1.
double F_func_1(double x);

/// Solution function for test case 1
double U_func_1(double x);

// *********************************************************************************************************************
// Test case 2 with f=0, g=1, h=1; solution is u(x) = 2-x
/// Force function for test case 2
double F_func_2(double x);

/// Solution function for test case 2
double U_func_2(double x);

// *********************************************************************************************************************
// Test case 3 with f=1, g=1, h=1; solution is u(x) = 2-x + 1/2(1-x^2)
/// Force function for test case 3
double F_func_3(double x);

/// Solution function for test case 3
double U_func_3(double x);

// *********************************************************************************************************************
// Test case 4; solution is u(x)  = sin(x) -->
// u''(x) = -sin(x), h=-u'(0)=-1, g=u(1)=sin(1)=0.841470984808
/// Force function for test case 4
double F_func_4(double x);

/// Solution function for test case 4
double U_func_4(double x);

// *********************************************************************************************************************
// Factory function to return a function pointer given name of a force function
// using FuncType = double (*FuncType) (double);
typedef double (*FuncType) (double);

// Initialize force table
void makeFuncTable();

/// Function returning a pointer to a force or solution function given its name
FuncType FuncByName(string funcName);
