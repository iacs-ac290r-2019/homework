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
#include "LinearSolve.hpp"

// *********************************************************************************************************************
// Declare the external LAPACK functions
extern "C" {
    // DGESV solves a general linear equation Ax = b by Gaussian elimination
    void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipivot,
                double *b, int *ldb, int *info);
}

// *********************************************************************************************************************
void lapack_fail() {
    cout << "LAPACK routine failed\n";
    exit(1);
}

// *********************************************************************************************************************
void solve_matrix(int n, double *A, double *x) {

    // Create the temporary memory that LAPACK needs
    int info, nrhs=1;
    int *ipiv = new int[n];

    // Make call to LAPACK routine
    dgesv_(&n, &nrhs, A, &n, ipiv, x, &n, &info);

    // Remove temporary memory
    delete [] ipiv;

    // Check for a non-zero value in info variable, indicating an error
    if (info!=0) 
    {
        lapack_fail();
    }
}
