#include <cstdio>

#include "lbm.hh"

int main() {
    char odir[]="lbm_mbc_Re001.out"; // output directory

    /** This is a simple Poiseuille flow.
     *  With pressure gradient on west and east.
     *  No-slip on north and south.
     *  initial condition: rho_0 = 1.0, ux_0 = Re*nu/D, uy_0 = 0.0
     *  Expected outcome: uy stays at 0.0, density stays at 1.0,
     *                    ux has parabolic behavior at vertical cut.
     *  No forcing. 
     *  Run 2000 time (relative) and output every 20 frame. */

    // Set up simulation parameters.
    double Re=0.01;      // Reynolds number
    double tau=1;    // Relaxation constant
    double D=60;       // Reference length    
    // Set up grid properties.
    int nx=200;        // Channel length
    int ny=60;         // Channel width 
    int w=30;          // Narrowing width
    int l=50;          // Narrowing length
    int flowtype=1;    // 0: No density input-output, 1: Inlet-Outlet
    int bctype=1;      // 0: Periodic, 1: Mixed
    double forcetype=0.;
    
    // Create the simulation domain with the specified parameters and dimensions
    lbm fl(Re,tau,D,nx,ny,odir);

    // Calculate the input velocity
    double u0=Re*fl.nu/D;

    // Call initialization functions to create simulation region and set up the initial condition
    fl.initialize(1.,u0,0.,flowtype);
    // double forcetype=0.0001;
    fl.solve(2000,100,bctype,forcetype,w,l);
    printf("Ran for 2000 iterations. Saved output every 20 frames. Finished!\n");
}