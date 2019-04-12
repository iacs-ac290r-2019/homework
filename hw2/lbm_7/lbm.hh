#ifndef LBM_HH
#define LBM_MM

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "lattice.hh"
#include "file.hh"

class lbm {
	public:
		/** Reynolds number. */
	    const double Re;
	    /** The relaxation time. */
	    const double tau;
	    /** kinematic viscosity. */
	    const double nu;
	    /** Characteristic length. */
	    const double D;
		/** The number of grid points in x direction, and padded with boundary buffer. */
		const int nx,nxb;
		/** The number of grid points in y direction, and padded with boundary buffer. */
		const int ny,nyb;		
		/** Pointers to lattice objects. */
		lattice *f;
		/** Boolean array of obstacles indicators. */
		bool *obs;

		/** Instantitate class object. */
		lbm(double Re_,double tau_,double D_,const int _nx,const int _ny,char *outdir_,const char *obsfile=NULL);
		/** Destory class object. */
		~lbm();

		/****************************** Initialization routine *******************************/
		/** Initialize routine called in main file. */
		void initialize(double macro_rho,double macro_ux,double macro_uy,int flowtype);
		/** Initialize hydro variables to all lattices. */
		void init_hydro(double macro_rho,double macro_ux,double macro_uy);
		/** Initialize flowtype. */
		void init_flow(int flowtype);
		/** Initialize Poiseuille flow. */
		void init_poiseuille();
		/** Initialize steady state Poiseuille flow. */
		void init_sstpoiseuille();
		/** Initialize populations. */
		void init_pop();

		/*********************************** Main routine ************************************/
		/** Main solve routine called in main file.*/
		void solve(int niters,int nout,int bctype,double forcetype,int w,int l);
		/** Set the boundary condition. */
		void bc(int bctype);
		/** Periodic boundary condition. */
		void pbc();
		/** Mixed boundary condition. */
		void mbc();
		/** Move the populations along eight directions. */
		void stream();
		/** Calculate hydro variables. */
		void hydro();
		/** Calculate the equilibrium population. */
		void equilibrium();
		/** Collide and calculate the population in the next timestep. */
		void collide();
		/** Add forcing. */
		void force(double forcetype);
		/** Add obstacle. */
		void obstacle(int i,int jbot,int jtop);
		/** Set up channel shape. */
		void set_channel(int w,int l);

		/***************************** Post-processing routine *******************************/
		/** Set output directory. */
		void set_dir(char *outdir_);
		/** Write the hydro quantities to file. */
		void output(int fr);
		/** Print f and f_eq. */
		void debug();
		/** Check on the sum of f. */
		void checkf();
	private:
		/** Output directory filename */
        char *outdir;
		/** Buffer for assembling output filenames */
        char *outbuf;
        /** Pointer to output file */
		FILE *fp;

};

#endif