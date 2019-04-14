#include "lbm.hh"

/** Constructor to the lbm class.
 * \param[in] nx the number of gridpoints in x direction.
 * \param[in] ny the number of gridpoints in y direction.
 * \param[in] outdir_ the output directory.
 * \param[in] obsfile an optional binary file of boolean flags denoting obstacles. */
lbm::lbm(double Re_,double tau_,double D_,const int nx_,const int ny_,char *outdir_,const char *obsfile)
	: Re(Re_), tau(tau_), nu((2*tau-1.)/6.), D(D_), 
	  nx(nx_), nxb(nx+2), ny(ny_), nyb(ny+2), 
	  f(new lattice[nxb*nyb]), obs(new bool[nxb*nyb]){
	// Set up output directory
	set_dir(outdir_);
	// Set up obstacle boolean array
    
}

/** Destory the class objects. */
lbm::~lbm() {
	delete [] obs;
	delete [] f;
	delete [] outdir;
}

/** Initialize hydro variables to all lattices. 
 * \param[in] macro_rho the macrosopic density.
 * \param[in] macro_ux the macrosopic velocity in the x direction.
 * \param[in] macro_uy the macrosopic velocity in the y direction. */
void lbm::init_hydro(double macro_rho,double macro_ux,double macro_uy) {
	for(int i=0;i<nxb;i++) {
		for(int j=0;j<nyb;j++) {
			int k=j*nxb+i;
			f[k].init_hydro(macro_rho,macro_ux,macro_uy);
		}
	}
}

/** Initialize flowtype.
 * \param[in] flowtype 1: Poiseuille flow. */
void lbm::init_flow(int flowtype) {
	if(flowtype==1) init_poiseuille();
	// Print out initial inputs
	printf("Poiseuille flow:\n");
	printf("Reynolds number: %.6f\n",Re);
	printf("Relaxation parameter: %.6f\n",tau);
	printf("Kinematic viscosity: %.6f\n",nu);
	printf("Characteristic length: %.0f\n",D);
	printf("Maximum steady velocity: %g\n",Re*nu/D);
	printf("Force: %g\n",8*nu*Re*nu/D/D/D);
	printf("Initial density at inlet: %g\n",f[nxb].rho);
	printf("Initial density at outlet: %g\n",f[nxb*nyb-1].rho);
}

/** Initialize Poiseuille flow. */
void lbm::init_poiseuille() {
	for(int j=1;j<=ny;j++) {
		f[j*nxb].rho=1+11./5000000/3600;
		// f[j*nxb+1].ux=Re*nu/D;
		f[(j+1)*nxb-1].rho=1-11./5000000/3600;
	}
}

/** Initialize populations. */
void lbm::init_pop() {
	for(int i=0;i<nxb;i++) {
		for(int j=0;j<nyb;j++) {
			int k=j*nxb+i;
			f[k].init_pop();
		}
	}
}

/** Set the boundary condition.
 * \param[in] bc_type the type of boundary condition. */
void lbm::bc(int bctype) {
	// Periodic boundary condition
	if(bctype==0) pbc();
	// Open boundary conidtion
	else if(bctype==1) obc();
}

/** Periodic boundary condition.
 *  West periodic, east periodic.
 *  North no-slip, south no-slip.  */
void lbm::pbc() {
	// West inlet
	for(int j=1;j<=ny;j++) {
		int k=j*nxb;
		f[k].f1=f[k+nx].f1;
		f[k].f5=f[k+nx].f5;
		f[k].f8=f[k+nx].f8;
	}
	// East outlet
	for(int j=1;j<=ny;j++) {
		int k=j*nxb+nx+1;
		f[k].f3=f[k-nx].f3;
		f[k].f6=f[k-nx].f6;
		f[k].f7=f[k-nx].f7;
	}
	// North solid
	for(int i=1;i<=nx;i++) {
		int k=nxb*(nyb-1)+i;
		f[k].f4=f[k-nxb].f2;
		f[k].f7=f[k-nxb-1].f5;
		f[k].f8=f[k-nxb+1].f6;
	}
	// South solid
	for(int i=1;i<=nx;i++) {
		int k=i;
		f[k].f2=f[k+nxb].f4;
		f[k].f5=f[k+nxb+1].f7;
		f[k].f6=f[k+nxb-1].f8;
	}
	// Northwest corner bounce-back
	f[nxb*(nyb-1)].f8=f[nxb*(nyb-1)-nxb+1].f6;
	// Northeast corner bounce-back
	f[nxb*nyb-1].f7=f[nxb*nyb-nxb-1-1].f5;
	// Southwest corner bounce-back
	f[0].f5=f[nxb+1].f7;
	// Southeast corner bounce-back
	f[nxb-1].f6=f[nxb+nxb-1-1].f8;
}

/** Open boundary condition.
 *  West inlet, east inlet.
 *  North no-slip, south no-slip. */
void lbm::obc() {
	// West inlet
	for(int j=1;j<=ny;j++) {
		int k=j*nxb;
		f[k].equilibrium();
		f[k].f1=f[k].feq1;
		f[k].f5=f[k].feq1;
		f[k].f8=f[k].feq1;
	}
	// East outlet
	// This version imposes pressure at the outlet
	// Hence no flux condition normal to the wall is not implemented
	for(int j=1;j<=ny;j++) {
		int k=j*nxb+nx+1;
		f[k].f3=f[k-1].f3;
		f[k].f6=f[k-1].f6;
		f[k].f7=f[k-1].f7;
	}
	// North solid
	for(int i=1;i<=nx;i++) {
		int k=nxb*(nyb-1)+i;
		f[k].f4=f[k-nxb].f2;
		f[k].f7=f[k-nxb-1].f5;
		f[k].f8=f[k-nxb+1].f6;
	}
	// South solid
	for(int i=1;i<=nx;i++) {
		int k=i;
		f[k].f2=f[k+nxb].f4;
		f[k].f5=f[k+nxb+1].f7;
		f[k].f6=f[k+nxb-1].f8;
	}
	// Northwest corner bounce-back
	f[nxb*(nyb-1)].f8=f[nxb*(nyb-1)-nxb+1].f6;
	// Northeast corner bounce-back
	f[nxb*nyb-1].f7=f[nxb*nyb-nxb-1-1].f5;
	// Southwest corner bounce-back
	f[0].f5=f[nxb+1].f7;
	// Southeast corner bounce-back
	f[nxb-1].f6=f[nxb+nxb-1-1].f8;

	init_poiseuille();
}


/** Stream the populations along eight directions. */
void lbm::stream() {
	lattice *z=new lattice[nxb*nyb];
    memcpy(z,f,nxb*nyb*sizeof(lattice));
	for(int i=1;i<=nx;i++) {
		for(int j=1;j<=ny;j++) {
			int k=j*nxb+i;
			f[k].f1=z[k-1].f1;
			f[k].f2=z[k-nxb].f2;
			f[k].f3=z[k+1].f3;
			f[k].f4=z[k+nxb].f4;
			f[k].f5=z[k-nxb-1].f5;
			f[k].f6=z[k-nxb+1].f6;
			f[k].f7=z[k+nxb+1].f7;
			f[k].f8=z[k+nxb-1].f8;
		}
	}
	delete [] z;
}

/** Calculate hydro variables. */
void lbm::hydro() {
	for(int i=1;i<=nx;i++) {
		for(int j=1;j<=ny;j++) {
			int k=j*nxb+i;
			f[k].hydro();
		}
	}
}

/** Calculate the equilibrium populations. */
void lbm::equilibrium() {
	for(int i=0;i<nxb;i++) {
		for(int j=0;j<nyb;j++) {
			int k=j*nxb+i;
			f[k].equilibrium();
		}
	}
}

/** Collide and calculate the population in the next timestep. */
void lbm::collide() {
	double omega=1./tau;
	for(int i=1;i<=nx;i++) {
		for(int j=1;j<=ny;j++) {
			int k=j*nxb+i;
			f[k].collide(omega);
		}
	}
}

/** Add forcing. */
void lbm::force(double forcetype) {
	for(int i=1;i<=nx;i++) {
		for(int j=1;j<=ny;j++) {
			int k=j*nxb+i;
			f[k].force(forcetype);
		}
	}
}

/** Add obstacle. 
 * \param[in] i the x position of the obstacle.
 * \param[in] jbot the bottom y position of the obstacle.
 * \param[in] jtop the top y position of the obstacle. */
void lbm::obstacle(int i,int jbot,int jtop) {
	for(int j=jbot;j<=jtop;j++) {
		int k=j*nxb+i;
		f[k].f1=f[k+1].f3;
		f[k].f5=f[k+nxb+1].f7;
		f[k].f8=f[k-nxb+1].f6;
		f[k].f3=f[k-1].f1;
		f[k].f7=f[k-nxb-1].f5;
		f[k].f6=f[k+nxb-1].f8;
	}
	f[jtop*nxb+i].f2=f[jtop*nxb+i+nxb].f4;
	f[jtop*nxb+i].f6=f[jtop*nxb+nxb-1].f8;
	f[jbot*nxb+i].f4=f[jbot*nxb+i-nxb].f2;
	f[jbot*nxb+i].f7=f[jbot*nxb+i-nxb-1].f5;
}

/** Set up the shape of the channel custom to the problem. 
 * \param[in] w the width of the narrowing.
 * \param[in] l the length of the narrowing. */
void lbm::set_channel(int w,int l) {
	int offset=(nx-l)/2;
	int h=(ny-w)/2;
	// Bottom two vertical
	obstacle(offset+1,1,h);
	obstacle(offset+l-1,1,h);
	// Top two vertical
	obstacle(offset+1,h+w,ny);
	obstacle(offset+l-1,h+w,ny);
	// North solid
	for(int i=offset+1;i<=offset+l-1;i++) {
		int k=nxb*(h+1)+i;
		f[k].f4=f[k-nxb].f2;
		f[k].f7=f[k-nxb-1].f5;
		f[k].f8=f[k-nxb+1].f6;
	}
	// South solid
	for(int i=offset+1;i<=offset+l-1;i++) {
		int k=nxb*(h+w-1)+i;
		f[k].f2=f[k+nxb].f4;
		f[k].f5=f[k+nxb+1].f7;
		f[k].f6=f[k+nxb-1].f8;
	}
}

/** Initialize routine called in main file.
 * \param[in] macro_rho the macrosopic density.
 * \param[in] macro_ux the macrosopic velocity in the x direction.
 * \param[in] macro_uy the macrosopic velocity in the y direction. */
void lbm::initialize(double macro_rho,double macro_ux,double macro_uy,int flowtype) {
	init_hydro(macro_rho,macro_ux,macro_uy);
	init_flow(flowtype);
	equilibrium();
	init_pop();
}

/** Main solver to step forward in time using the lattice Boltzmann method.
 * \param[in] niters the number of iterations.
 * \param[in] nout the number of equally-spaced output frames.
 * \param[in] bc_type the type of boundary condition.
 * \param[in] forcetype the type and value of external force.
 * \param[in] w the width of the narrowing.
 * \param[in] l the length of the narrowing. */
void lbm::solve(int niters,int nout,int bctype,double forcetype,int w,int l) {
	int skip=niters/(nout-1);
	int fr=0;
	output(fr);
	fr++;
	for(int t=1;t<niters;t++) {
		// LBM routine
		bc(bctype);
		set_channel(w,l);
	    stream();
	    hydro();
	    equilibrium(); 
	    collide();
	    force(forcetype);	 
		// Output to files
		if(t%skip==0) {
			output(fr);
			fr++;
		}
	}
}

/** Set up output directory 
 * \param[in] outdir_ the output directory name */
void lbm::set_dir(char *outdir_) {
	size_t l=strlen(outdir_)+1;
	outdir=new char[2*l+32];
    memcpy(outdir,outdir_,sizeof(char)*l);
    outbuf=outdir+l;
    // Create output directory, if it doesn't exist.
	mkdir(outdir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
}

/** Output the hydro quantities to files
 * \param[in] fr the frame number */
void lbm::output(int fr) {
	// Header line prints the timestep and system dimensions
	sprintf(outbuf,"%s/fr_%04d.txt",outdir,fr);
	fp=safe_fopen(outbuf,"w");
	fprintf(fp,"0 0 0 %d %d %d\n",fr,nx,ny);
	// Output the density, x-velocity, and y-velocity
	for(int i=1;i<=ny;i++) {
		for(int j=1;j<=nx;j++) {
			fprintf(fp,"%d %d %g %g %g %g\n",i,j,f[i*nxb+j].rho,f[i*nxb+j].ux,f[i*nxb+j].uy, f[i*nxb+j].p);
		}
	}
	fclose(fp);
}

void lbm::debug() {
	for(int i=0;i<nyb;i++) {
		for(int j=0;j<nxb;j++) {
			printf("\nDebug Lattice %d\n",i*nxb+j);
			printf("Debug f0: %g, feq0: %g\n",f[i*nxb+j].f0,f[i*nxb+j].feq0);
			printf("Debug f1: %g, feq1: %g\n",f[i*nxb+j].f1,f[i*nxb+j].feq1);
			printf("Debug f2: %g, feq2: %g\n",f[i*nxb+j].f2,f[i*nxb+j].feq2);
			printf("Debug f3: %g, feq3: %g\n",f[i*nxb+j].f3,f[i*nxb+j].feq3);
			printf("Debug f4: %g, feq4: %g\n",f[i*nxb+j].f4,f[i*nxb+j].feq4);
			printf("Debug f5: %g, feq5: %g\n",f[i*nxb+j].f5,f[i*nxb+j].feq5);
			printf("Debug f6: %g, feq6: %g\n",f[i*nxb+j].f6,f[i*nxb+j].feq6);
			printf("Debug f7: %g, feq7: %g\n",f[i*nxb+j].f7,f[i*nxb+j].feq7);
			printf("Debug f8: %g, feq8: %g\n",f[i*nxb+j].f8,f[i*nxb+j].feq8);
		}
	}
}

void lbm::checkf() {
	double sumxf=0;
	double sumxfeq=0;
	for(int i=1;i<=ny;i++) {
		for(int j=1;j<=nx;j++) {
			double sumf=0;
			double sumfeq=0;
			sumf=f[i*nxb+j].f0+f[i*nxb+j].f1+f[i*nxb+j].f2+f[i*nxb+j].f3+f[i*nxb+j].f4
				+f[i*nxb+j].f5+f[i*nxb+j].f6+f[i*nxb+j].f7+f[i*nxb+j].f8;
			sumfeq=f[i*nxb+j].feq0+f[i*nxb+j].feq1+f[i*nxb+j].feq2+f[i*nxb+j].feq3+f[i*nxb+j].feq4
				+f[i*nxb+j].feq5+f[i*nxb+j].feq6+f[i*nxb+j].feq7+f[i*nxb+j].feq8;
			printf("Lattice %d, sum of f: %g, sum of f_eq: %g\n",i*nxb+j,sumf,sumfeq);
			sumxf+=sumf;
			sumxfeq+=sumfeq;
		}
	}
	printf("Sum x of sum of f: %g, sum x of sum of f_eq: %g\n",sumxf,sumxfeq);
}