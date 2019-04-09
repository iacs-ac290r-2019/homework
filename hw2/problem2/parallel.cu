/* Daniel Willen, 2019
 *
 * Solve the transient heat conduction problem with homogeneous Dirichlet
 *  boundary conditions:
 *
 *    u(x={0,L}) = u(y={0,L}) = 0
 *
 *  and initial condition:
 *
 *    u(x,y,0) = sin(x) * sin(y)
 *
 *  on the domain 0 <= x,y <= L, with L = pi.
 *
 * This program solves the above problem on a single GPU with the Jacobi method.
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>

#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

#define PI 3.14159265358979323846
#define MAX_THREADS_DIM 16        // Note that this depends on the hardware

/* Note on the structure of this file:
 *  - Cuda device constant memory declarations are at the top
 *  - Functions definitions are in the middle. Functions include:
 *  - - parse_cmdline: Read command-line arguments for domain size
 *  - - jacobi_solver: Advance the soln to the next time step using Jacobi
 *  - - check_error:   Calculate the error b/t the numeric and analytic solns
 *  - The `main' function is at the bottom
 *
 *  Note that it is good practice to use header files and break functions out
 *   into separate files. This has not been done here for simplicity.
 */

/*** Auxiliary Functions ***/

/* Read the command line inputs */
// - argv[0] is the program name
// - argv[1] is the first input (number of points)
int parse_cmdline(int argc, char *argv[]) {
  int nx;
  if (argc == 2) {
    nx = atoi(argv[1]); // Number of grid points
		if (nx < MAX_THREADS_DIM) {
			printf("Expecting a number of grid cells in one dimension to be at least %d\n", MAX_THREADS_DIM);
			exit(EXIT_FAILURE);
		}

    printf("Grid is %d by %d\n\n", nx, nx);
  } else {
    printf("Input error. Run like: \n\n");
    printf("  $ ./parallel.c n\n\n");
    printf("  where n is the number of grid cells in one dimension\n");
    exit(EXIT_FAILURE);
  }
  return nx;
}

/*******************************************************************************
 * Step IV: Launch the GPU kernel to advance to the next time step with the    *
 *          Jacobi method here.                                                *
 ******************************************************************************/
__global__ void computeNextJacobiStep(int nx, int ny, double pref, double* _u, double* _u_new) {
	int ti = blockDim.x * blockIdx.x + threadIdx.x;
	int tj = blockDim.y * blockIdx.y + threadIdx.y;

	if (ti < (nx-1) && ti > 0 && tj < (ny-1) && tj > 0) {
		double leftTerm = _u[tj*nx + ti];
		double rightTerm = pref * (
			_u[tj*nx + (ti+1)] +
			_u[tj*nx + (ti-1)] +
			_u[(tj+1)*nx + ti] +
			_u[(tj-1)*nx + ti] -
			4*_u[tj*nx + ti]
		);
		
		_u_new[tj*nx + ti] = leftTerm + rightTerm;
	}
}

/******************************************************************************
 * Step V: Launch the GPU kernel to calculate the error at each grid point    *
 *         here.                                                              *
 *****************************************************************************/
__global__ void computeJacobiError(int nx, int ny, double D, double t, double* _u, double* _error) {
	int ti = blockDim.x * blockIdx.x + threadIdx.x;
	int tj = blockDim.y * blockIdx.y + threadIdx.y;

	if (ti < (nx-1) && ti > 0 && tj < (ny-1) && tj > 0) {
		double discretizedValue = _u[tj*nx + ti];
		double analyticalValue = sin(ti)*sin(tj)*exp(-2*D*t);

		_error[tj*nx + ti] = abs(discretizedValue - analyticalValue);
	}
}

/*** Main Function ***/
int main(int argc, char *argv[])
{
  /* Variable declaration */
  double Lx = PI;           // Domain length in x-direction
  double Ly = PI;           // Domain length in y-direction
  double D = 1.;            // Diffusion constant

  int nx, ny;               // Grid points (grid cells + 1)
  double dx, dy;            // Grid spacing
  double dt;                // Time step size
  double sim_time;          // Length of sim time, arbitrary for simplicity
  double pref;              // Pre-factor in the Jacobi method

  double error = 0.;        // Mean percent-difference at each grid point
  error = error;            // To prevent compiler warning

  /* Parse command-line for problem size */
  nx = parse_cmdline(argc, argv);
  ny = nx;                  // Assume a square grid

  /* Initialize variables */
  dx = Lx / (nx - 1);       // Cell width in x-direction
  dy = Ly / (ny - 1);       // Cell width in y-direction
  dt = 0.25*dx*dy/D;        // Limited by diffusive stability
  sim_time = 0.5*Lx*Ly/D;   // Arbitrary simulation length
  pref = D*dt/(dx*dx);      // Jacobi pre-factor

  /*****************************************************************************
   * Step I: Declare, allocate, and initialize memory for the field variable   *
   *         u on the CPU.                                                     *
   ****************************************************************************/
	double* u = (double*) malloc(nx*ny * sizeof(double));
	for (int j = 0; j < ny; ++j) {
		for (int i = 0; i < nx; ++i) {
			u[j*nx + i] = sin(i) * sin(j); 
		}
	}

  /*****************************************************************************
   * Step II: Declare and allocate GPU memory for _u, _u_new, and _error. Copy *
   *          the initial condition to the GPU.                                *
   ****************************************************************************/
	double *_u, *_u_new, *_error;
	cudaMalloc(&_u, nx*ny * sizeof(double));
	cudaMemcpy(_u, u, nx*ny * sizeof(double), cudaMemcpyHostToDevice);
	cudaMalloc(&_u_new, nx*ny * sizeof(double));
	cudaMalloc(&_error, nx*ny * sizeof(double));

  // Set the new soln and error to 0
  cudaMemset(_u_new, 0., nx*ny * sizeof(double));
  cudaMemset(_error, 0., nx*ny * sizeof(double));

  // Create thrust pointers to device memory for error calculation
  thrust::device_ptr<double> t_error(_error);

  /*****************************************************************************
   * Step III: Set up the kernel execution configuration for the domain based  *
   *           on the input domain size and the MAX_THREADS_DIM variable.      *
   ****************************************************************************/
	int tx = MAX_THREADS_DIM;
	int ty = MAX_THREADS_DIM;

	int bx = (int) ceil((double) nx / tx);
	int by = (int) ceil((double) ny / ty);

	dim3 dimBlocks(tx, ty);
	dim3 numBlocks(bx, by);

  /***************************/
  /* Main Time-Stepping Loop */
  /***************************/
  for (double time = 0.; time <= sim_time; time += dt) {
    /***************************************************************************
     * Step IV: Launch the GPU kernel to advance to the next time step with    *
     *          the Jacobi method here.                                        *
     **************************************************************************/
		computeNextJacobiStep<<<numBlocks, dimBlocks>>>(nx, ny, pref, _u, _u_new);
		cudaDeviceSynchronize();

    /***************************************************************************
     * Step V: Launch the GPU kernel to calculate the error at each grid point *
     *         here.                                                           *
     **************************************************************************/
		computeJacobiError<<<numBlocks, dimBlocks>>>(nx, ny, D, time, _u, _error);
		cudaDeviceSynchronize();

    // Use thrust to do a parallel reduction on the error
    error = thrust::reduce(t_error, t_error + nx*ny, 0., thrust::plus<double>());
    printf("Error at t* = %.5lf is %e\n", time*D/(Lx*Lx), error/(nx*ny));

    // Copy new soln to old. This also blocks to ensure computations are finished.
    cudaMemcpy(_u, _u_new, nx*ny * sizeof(double), cudaMemcpyDeviceToDevice);
  }

  /*****************************************************************************
   * Step VI: Copy the memory back to the CPU.                                 *
   ****************************************************************************/
	cudaMemcpy(u, _u, nx*ny * sizeof(double), cudaMemcpyDeviceToHost);

  /*****************************************************************************
   * Step I and Step II: Free the memory that you declared and allocated       *
   *                     earlier in the program.                               *
   ****************************************************************************/
	free(u);
	cudaFree(_u);
	cudaFree(_u_new);
	cudaFree(_error);

  return EXIT_SUCCESS;
}

