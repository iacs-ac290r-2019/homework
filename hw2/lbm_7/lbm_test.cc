#include <cstdio>

#include "lbm.hh"

int main() {
    char odir[]="cyl_Re80.out"; // output directory
    char buf[]="cyl24.bin";  // file encoding barrier flags

    /** This is a simple Poiseuille flow.
     *  Periodic boundary condition on west and east.
     *  No-slip on north and south.
     *  initial condition: rho_0 = 1.0, ux_0 = 0.0, uy_0 = 0.0
     *  Expected outcome: uy stays at 0.0, density stays at 1.0,
     *                    ux has parabolic behavior at vertical cut.
     *  forcing gives the velocity. 
     *  Run 1000 time (relative) and output every 100 frame. */

    // Set up simulation parameters.
    double Re=0.01;      // Reynolds number
    double tau=1;    // Relaxation constant
    double D=1;	   // Reference length
    // Set up grid properties.
    int nx=200;		   // Channel length
    int ny=40;		   // Channel width	
    int flowtype=0;    // 0: Static flow, 1: Poiseuille flow
    int bctype=0;      // 0: Periodic, 1: Inlet-Outlet
    
    // Create the simulation domain with the specified parameters and dimensions
    lbm fl(Re,tau,D,nx,ny,odir);

    // double forcetype=1./8*0.1*fl.nu/nx/nx; // Poiseuille flow, NOT USED THOUGH

    // Call initialization functions to create simulation region and set up the initial condition
    fl.initialize(1.,0.,0.,flowtype);
    double forcetype=0.0001;
    fl.solve(1000,100,bctype,forcetype);
    printf("Great!\n");
}
// 