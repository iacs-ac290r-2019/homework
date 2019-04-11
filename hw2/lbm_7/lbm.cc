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
			f[k].rho=macro_rho;
			f[k].ux=macro_ux;
			f[k].uy=macro_uy;
		}
	}
}

/** Initialize flowtype.
 * \param[in] flowtype 0: Poiseuille flow. */
void lbm::init_flow(int flowtype) {
	if(flowtype==1) init_poiseuille();
}

/** Initialize Poiseuille flow. */
void lbm::init_poiseuille() {
	// Used later for a physically accurate Poiseuille flow
	// double Ly=static_cast<double>(ny),Lx=static_cast<double>(nx);
	// double dp=36*Re*nu*nu*Lx/D/Ly/Ly;
	for(int i=0;i<nyb;i++) {
		f[i*nxb+1].rho=1+0.01;
		f[(i+1)*nxb-2].rho=1-0.01;
	}
	// Print out initial inputs
	printf("Poiseuille flow:\n");
	printf("Reynolds number: %.6f\n",Re);
	printf("Relaxation parameter: %.6f\n",tau);
	printf("Kinematic viscosity: %.6f\n",nu);
	printf("Characteristic length: %.0f\n",D);
	printf("Initial density at inlet: %.6f\n",f[nxb+1].rho);
	printf("Initial density at outlet: %.6f\n",f[nxb*nyb-2].rho);
	printf("Centerline velocity: %.6f\n",3*Re*nu/2/D);
}

/** Initialize steady state Poiseuille flow. */
// void lbm::init_sstpoiseuille() {
// 	double u,y,p,Ly=static_cast<double>(ny-1),Lx=static_cast<double>(nx-1);
// 	double dp=36*Re*nu*nu*Lx/D/Ly/Ly;
// 	// East and west
// 	for(int i=0;i<nyb;i++) {
// 		y=1-1./Ly*i;
// 		u=6*Re*nu/D*(y-y*y);
// 		f[i*nxb+1].rho=1+dp/2;
// 		f[i*nxb+1].ux=u;
// 		f[(i+1)*nxb-2].rho=1-dp/2;
// 		f[(i+1)*nxb-2].ux=u;
// 	}
// 	// Top and bottom
// 	for(int j=1;j<=nx;j++) {
// 		p=1+dp*(1./2-j/Lx);
// 		f[j].uy=0;
// 		f[j].rho=p;
// 		f[nxb*(nyb-1)+j].uy=0;
// 		f[nxb*(nyb-1)+j].rho=p;
// 	}
// }

/** Initialize populations. */
void lbm::init_pop() {
	lattice *z=new lattice[nxb*nyb];
    memcpy(z,f,nxb*nyb*sizeof(lattice));
	for(int i=0;i<nxb;i++) {
		for(int j=0;j<nyb;j++) {
			int k=j*nxb+i;
			f[k].f0=z[k].feq0;
			f[k].f1=z[k].feq1;
			f[k].f2=z[k].feq2;
			f[k].f3=z[k].feq3;
			f[k].f4=z[k].feq4;
			f[k].f5=z[k].feq5;
			f[k].f6=z[k].feq6;
			f[k].f7=z[k].feq7;
			f[k].f8=z[k].feq8;
		}
	}
	delete [] z;
}

/** Set the boundary condition.
 * \param[in] bc_type the type of boundary condition. */
void lbm::bc(int bctype) {
	// Periodic boundary condition
	if(bctype==0) pbc();
	// Mixed boundary conidtion
	else if(bctype==1) mbc();
}

/** Periodic boundary condition.
 *  West periodic, east periodic.
 *  North no-slip, south no-slip.  */
void lbm::pbc() {
	lattice *z=new lattice[nxb*nyb];
    memcpy(z,f,nxb*nyb*sizeof(lattice));
	// West inlet
	for(int j=1;j<=ny;j++) {
		int k=j*nxb;
		f[k].f1=z[k+nx].f1;
		f[k].f5=z[k+nx].f5;
		f[k].f8=z[k+nx].f8;
		// printf("West inlet done\n");
	}
	// East outlet
	for(int j=1;j<=ny;j++) {
		int k=j*nxb+nx+1;
		f[k].f3=z[k-nx].f3;
		f[k].f6=z[k-nx].f6;
		f[k].f7=z[k-nx].f7;
		// printf("East outlet done\n");
	}
	// North solid
	for(int i=1;i<=nx;i++) {
		int k=nxb*(nyb-1)+i;
		f[k].f4=z[k-nxb].f2;
		f[k].f7=z[k-nxb-1].f5;
		f[k].f8=z[k-nxb+1].f6;
	}
	// South solid
	for(int i=1;i<=nx;i++) {
		int k=i;
		f[k].f2=z[k+nxb].f4;
		f[k].f5=z[k+nxb+1].f7;
		f[k].f6=z[k+nxb-1].f8;
	}
	// Northwest corner bounce-back
	f[nxb*(nyb-1)].f8=z[nxb*(nyb-1)-nxb+1].f6;
	// Northeast corner bounce-back
	f[nxb*nyb-1].f7=z[nxb*nyb-nxb-1-1].f5;
	// Southwest corner bounce-back
	f[0].f5=z[nxb+1].f7;
	// Southeast corner bounce-back
	f[nxb-1].f6=z[nxb+nxb-1-1].f8;
	delete [] z;
}

/** Mixed boundary condition.
 *  West inlet, east inlet.
 *  North no-slip, south no-slip.  */
void lbm::mbc() {
	// West inlet
	for(int j=1;j<=ny;j++) {
		int k=j*nxb;
		f[k].f1=f[k+1].f1;
		f[k].f5=f[k+1].f5;
		f[k].f8=f[k+1].f8;
		// printf("West inlet done\n");
	}
	// East outlet
	for(int j=1;j<=ny;j++) {
		int k=j*nxb+nx+1;
		f[k].f3=f[k-1].f3;
		f[k].f6=f[k-1].f6;
		f[k].f7=f[k-1].f7;
		// printf("East outlet done\n");
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
			f[k].rho=f[k].f0+f[k].f1+f[k].f2+f[k].f3+f[k].f4+f[k].f5+f[k].f6+f[k].f7+f[k].f8;
			f[k].ux=(f[k].f1-f[k].f3+f[k].f5-f[k].f6-f[k].f7+f[k].f8)/f[k].rho;
			f[k].uy=(f[k].f2-f[k].f4+f[k].f5+f[k].f6-f[k].f7-f[k].f8)/f[k].rho;
			// printf("Lattice %d, rho: %g, ux: %g, uy: %g\n",i*nxb+j,f[i*nxb+j].rho,f[i*nxb+j].ux,f[i*nxb+j].uy);
		}
	}
}

/** Calculate the equilibrium populations. */
void lbm::equilibrium() {
	// Sound speed squared
	// const double cs2=1./3,cs22=2*cs2,cssq=2./9;
	const double cs2=1./3,cs4=1./9;
	// Weights for D2Q9
	const double w0=4./9,w1=1./9,w2=1./36;
	// Traverse the simulation region with buffer to calculate f_eq
	lattice *z=new lattice[nxb*nyb];
    memcpy(z,f,nxb*nyb*sizeof(lattice));
	for(int i=0;i<nxb;i++) {
		for(int j=0;j<nyb;j++) {
			int k=j*nxb+i;
			double rho=z[k].rho;
			double u=z[k].ux;
			double v=z[k].uy;
			double ui=u/cs2;
			double vi=v/cs2;
			double u2=u*u/2/cs4;
			double v2=v*v/2/cs4;
			double sumsq=(u*u+v*v)/2/cs2;
			double sumsq2=sumsq*(1.-cs2)/cs2;
			double uv=ui*vi;
			f[k].feq0=rho*w0*(1.-sumsq);
			f[k].feq1=rho*w1*(1.-sumsq+ui+u2);
			f[k].feq2=rho*w1*(1.-sumsq+vi+v2);
			f[k].feq3=rho*w1*(1.-sumsq-ui+u2);
			f[k].feq4=rho*w1*(1.-sumsq-vi+v2);
			f[k].feq5=rho*w2*(1.+sumsq2+ui+vi+uv);
			f[k].feq6=rho*w2*(1.+sumsq2-ui+vi-uv);
			f[k].feq7=rho*w2*(1.+sumsq2-ui-vi+uv);
			f[k].feq8=rho*w2*(1.+sumsq2+ui-vi-uv);
		}
	}
	delete [] z;
}

/** Collide and calculate the population in the next timestep. */
void lbm::collide() {
	double omega=1./tau;
	lattice *z=new lattice[nxb*nyb];
    memcpy(z,f,nxb*nyb*sizeof(lattice));
	for(int i=1;i<=nx;i++) {
		for(int j=1;j<=ny;j++) {
			int k=j*nxb+i;
			f[k].f0=(1.-omega)*z[k].f0+omega*z[k].feq0;
			f[k].f1=(1.-omega)*z[k].f1+omega*z[k].feq1;
			f[k].f2=(1.-omega)*z[k].f2+omega*z[k].feq2;
			f[k].f3=(1.-omega)*z[k].f3+omega*z[k].feq3;
			f[k].f4=(1.-omega)*z[k].f4+omega*z[k].feq4;
			f[k].f5=(1.-omega)*z[k].f5+omega*z[k].feq5;
			f[k].f6=(1.-omega)*z[k].f6+omega*z[k].feq6;
			f[k].f7=(1.-omega)*z[k].f7+omega*z[k].feq7;
			f[k].f8=(1.-omega)*z[k].f8+omega*z[k].feq8;
		}
	}
	delete [] z;
}

/** Add forcing. */
void lbm::force(double forcetype) {
	// Sound speed squared
	const double cs2=1./3;
	// Weights for D2Q9
	const double w1=1./9,w2=1./36;
	lattice *z=new lattice[nxb*nyb];
    memcpy(z,f,nxb*nyb*sizeof(lattice));
	for(int i=1;i<=nx;i++) {
		for(int j=1;j<=ny;j++) {
			int k=j*nxb+i;
			double rho=f[k].rho;
			f[k].f1=z[k].f1+w1*forcetype/cs2*rho;
			f[k].f5=z[k].f5+w2*forcetype/cs2*rho;
			f[k].f8=z[k].f8+w2*forcetype/cs2*rho;
			f[k].f3=z[k].f3-w1*forcetype/cs2*rho;
			f[k].f6=z[k].f6-w2*forcetype/cs2*rho;
			f[k].f7=z[k].f7-w2*forcetype/cs2*rho;
			// f[k].f1=f[k].f1+forcetype;
			// f[k].f5=f[k].f5+forcetype;
			// f[k].f8=f[k].f8+forcetype;
			// f[k].f3=f[k].f3-forcetype;
			// f[k].f6=f[k].f6-forcetype;
			// f[k].f7=f[k].f7-forcetype;
		}
	}
	delete [] z;
}

/** Add obstacle. */
void lbm::obstacle() {
	int nobst=10;
	int j=nx/4;
	int ibot=ny/2-nobst/2;
	int itop=ny/2+nobst/2+1;
	for(int i=ibot;i<=itop;i++) {
		f[i*nxb+j].f1=f[i*nxb+j+1].f3;
		f[i*nxb+j].f5=f[(i+1)*nxb+j+1].f7;
		f[i*nxb+j].f8=f[(i+1)*nxb+j-1].f6;
		f[i*nxb+j].f3=f[i*nxb+j-1].f1;
		f[i*nxb+j].f7=f[(i-1)*nxb+j].f5;
		f[i*nxb+j].f6=f[(i-1)*nxb+j].f8;
	}
	f[itop*nxb+j].f2=f[itop*nxb+j+1].f4;
	f[itop*nxb+j].f6=f[(itop-1)*nxb+j+1].f8;
	f[ibot*nxb+j].f4=f[ibot*nxb+j-1].f2;
	f[ibot*nxb+j].f7=f[(ibot-1)*nxb+j-1].f5;
}

/** Initialize routine called in main file.
 * \param[in] macro_rho the macrosopic density.
 * \param[in] macro_ux the macrosopic velocity in the x direction.
 * \param[in] macro_uy the macrosopic velocity in the y direction. */
void lbm::initialize(double macro_rho,double macro_ux,double macro_uy,int flowtype) {
	init_hydro(macro_rho,macro_ux,macro_uy);
	equilibrium();
	init_pop();
}

/** Main solver to step forward in time using the lattice Boltzmann method.
 * \param[in] niters the number of iterations.
 * \param[in] nout the number of equally-spaced output frames.
 * \param[in] bc_type the type of boundary condition.
 * \param[in] forcetype the type and value of external force. */
void lbm::solve(int niters,int nout,int bctype,double forcetype) {
	int skip=niters/(nout-1);
	int fr=0;
	// Output the initial frame
	output(fr);
	fr++;
	for(int t=1;t<niters;t++) {
		// printf("\n\n\n\nIteration %d\n",t);
		// LBM routine
		// init_poiseuille();
		bc(bctype);
		// printf("\nIteration %d Check f after bc()\n",t);
		// checkf();
	    stream();
	    // printf("\nIteration %d Check f after stream()\n",t);
		// checkf();
	    hydro();
	    // printf("\nIteration %d Check f after hydro()\n",t);
		// checkf();
	    equilibrium(); 
	    // printf("\nIteration %d Check f after equilibrium()\n",t);
		// checkf();
	    collide();
	    // printf("\nIteration %d Check f after collide()\n",t);
		// checkf();
	    force(forcetype);
	    // printf("\nIteration %d Check f after force()\n",t);
		// checkf();
	    // obstacle();
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
	fprintf(fp,"0 0 %d %d %d\n",fr,nx,ny);
	// Output the density, x-velocity, and y-velocity
	for(int i=1;i<=ny;i++) {
		for(int j=1;j<=nx;j++) {
			fprintf(fp,"%d %d %g %g %g\n",i,j,f[i*nxb+j].rho,f[i*nxb+j].ux,f[i*nxb+j].uy);
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