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
		/** Destory function. */
		virtual ~lattice() {}

		/** Initialize hydro variables. 
		 * \param[in] macro_rho the macrosopic density.
		 * \param[in] macro_ux the macrosopic velocity in the x direction.
		 * \param[in] macro_uy the macrosopic velocity in the y direction. */
		inline void init_hydro(double macro_rho,double macro_ux,double macro_uy) {
			rho=macro_rho;
			ux=macro_ux;
			uy=macro_uy;
		}

		/** Initialize populations with f_eq. */
		inline void init_pop() {
			f0=feq0; f1=feq1; f2=feq2; f3=feq3; f4=feq4;
			f5=feq5; f6=feq6; f7=feq7; f8=feq8;
		}

		/** Calculate hydro variables rho, ux and uy. */
		inline void hydro() {
			rho=f0+f1+f2+f3+f4+f5+f6+f7+f8;
			ux=(f1+f5+f8-f3-f6-f7)*1./rho;
			uy=(f2+f5+f6-f4-f7-f8)*1./rho;
		}

		/** Calculate the equilibrium populations. */
		inline void equilibrium() {
			// Sound speed constants
			const double cs2=1./3,cs4=1./9;
			// Weights for D2Q9
			const double w0=4./9,w1=1./9,w2=1./36;
			double ui=ux/cs2;
			double vi=uy/cs2;
			double u2=ux*ux/2/cs4;
			double v2=uy*uy/2/cs4;
			double sumsq=(ux*ux+uy*uy)/2/cs2;
			double sumsq2=sumsq*(1.-cs2)/cs2;
			double uv=ui*vi;
			feq0=rho*w0*(1.-sumsq);
			feq1=rho*w1*(1.-sumsq+ui+u2);
			feq2=rho*w1*(1.-sumsq+vi+v2);
			feq3=rho*w1*(1.-sumsq-ui+u2);
			feq4=rho*w1*(1.-sumsq-vi+v2);
			feq5=rho*w2*(1.+sumsq2+ui+vi+uv);
			feq6=rho*w2*(1.+sumsq2-ui+vi-uv);
			feq7=rho*w2*(1.+sumsq2-ui-vi+uv);
			feq8=rho*w2*(1.+sumsq2+ui-vi-uv);
		}

		/** Collide and calculate the population in the next timestep. 
		 * \param[in] omega 1/tau. */
		inline void collide(double omega) {
			f0=(1.-omega)*f0+omega*feq0;
			f1=(1.-omega)*f1+omega*feq1;
			f2=(1.-omega)*f2+omega*feq2;
			f3=(1.-omega)*f3+omega*feq3;
			f4=(1.-omega)*f4+omega*feq4;
			f5=(1.-omega)*f5+omega*feq5;
			f6=(1.-omega)*f6+omega*feq6;
			f7=(1.-omega)*f7+omega*feq7;
			f8=(1.-omega)*f8+omega*feq8;
		}

		/** Add forcing. */
		inline void force(double forcetype) {
			// Sound speed squared
			const double cs2=1./3;
			// Weights for D2Q9
			const double w1=1./9,w2=1./36;
			f1=f1+w1*forcetype/cs2*rho;
			f5=f5+w2*forcetype/cs2*rho;
			f8=f8+w2*forcetype/cs2*rho;
			f3=f3-w1*forcetype/cs2*rho;
			f6=f6-w2*forcetype/cs2*rho;
			f7=f7-w2*forcetype/cs2*rho;
		}
};

#endif