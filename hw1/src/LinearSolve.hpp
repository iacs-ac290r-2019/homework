/**
 *  @file LinearSolve.hpp
 *  @brief Solve linear equations by wrapping calls to LAPACK routines.
 * 
 *  SOURCE ATTRIBUTION
 *  Applied Math 225 Examples 2a includes a superset of the wrappers used here.
 *  This is NOT original work!
 * 
 *  @author Michael S. Emanuel
 *  @date 2019-02-22
 */

// *********************************************************************************************************************
#pragma once

// *********************************************************************************************************************
// Libraries
#include <iostream>
    using std::cout;

// *********************************************************************************************************************
/** Prints a generic error message in cases when LAPACK reports an error. */
void lapack_fail();

/** Solves the matrix system Ax=b.
 * \param[in] n the dimension of the matrix.
 * \param[in] A the matrix terms.
 * \param[in] x the source vector, which will be replaced by the solution when the routine exits. */
void solve_matrix(int n, double *A, double *x);
