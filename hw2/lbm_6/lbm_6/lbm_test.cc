#include <cstdio>

#include "lbm.hh"

int main() {
    char odir[]="cyl_Re80.out"; // output directory
    char buf[]="cyl24.bin";  // file encoding barrier flags

    // Set up simulation parameters.
    double Re=0.01;      // Reynolds number
    double tau=0.6;    // Relaxation constant
    double D=1;	   // Reference length
    // Set up grid properties.
    int nx=200;		   // Channel length
    int ny=40;		   // Channel width	
    int flowtype=0;    // 0: Poiseuille flow
    int bctype=0;      // 0: Mixed, else: to be implemented
    
    // Create the simulation domain with the specified parameters and dimensions
    lbm fl(Re,tau,D,nx,ny,odir);

    double forcetype=1./8*0.1*fl.nu/nx/nx; // Poiseuille flow, NOT USED THOUGH
    printf("Force value: %g\n",forcetype);

    // Call initialization functions to create simulation region and set up the initial condition
    double u=tau*Re/nx;
    printf("U value: %g\n",u);
    u=0.;
    fl.initialize(1.,u,0.,flowtype);
    
    // Main routine
    forcetype=0.;
    fl.solve(100,100,bctype,forcetype);
    printf("Great!\n");
}
