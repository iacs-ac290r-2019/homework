#ifndef LATTICE_HH
#define LATTICE_HH

#include <cstdio>
#include <cmath>

class lattice {
	public:
		/** Hydro variables. */
		double rho,ux,uy;
		/** Population at current state, equilibrium at the lattice */
		double f0,feq0;
		/** Population at current state, equilibrium streaming east */
		double f1,feq1;
		/** Population at current state, equilibrium streaming north */
		double f2,feq2;
		/** Population at current state, equilibrium streaming west */
		double f3,feq3;
		/** Population at current state, equilibrium streaming south */
		double f4,feq4;
		/** Population at current state, equilibrium streaming northeast */
		double f5,feq5;
		/** Population at current state, equilibrium streaming northwest */
		double f6,feq6;
		/** Population at current state, equilibrium streaming southwest */
		double f7,feq7;
		/** Population at current state, equilibrium streaming southeast */
		double f8,feq8;
		/** Boundary flag. */
		bool bf;
};

#endif