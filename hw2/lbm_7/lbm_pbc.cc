#include <cstdio>

#include "lbm.hh"

int main() {
    char odir[]="lbm_pbc_Re001.out"; // output directory

    /** This is a simple Poiseuille flow.
     *  Periodic boundary condition on west and east.
     *  No-slip on north and south.
     *  initial condition: rho_0 = 1.0, ux_0 = 0.0, uy_0 = 0.0
     *  Expected outcome: uy stays at 0.0, density stays at 1.0,
     *                    ux has parabolic behavior at vertical cut.
     *  forcing gives the velocity. 
     *  Run 15000 time (relative) and output every 150 frame. */

    // Set up simulation parameters.
    double Re=0.01;      // Reynolds number
    double tau=1;    // Relaxation constant
    double D=60;	   // Reference length
    // Set up grid properties.
    int nx=200;		   // Channel length
    int ny=60;		   // Channel width	
    int w=30;          // Narrowing width
    int l=50;          // Narrowing length
    int flowtype=0;    // 0: Static flow, 1: Poiseuille flow
    int bctype=0;      // 0: Periodic, 1: Inlet-Outlet
    double forcetype=0.0001; // Poiseuille flow, NOT USED THOUGH
    
    // Create the simulation domain with the specified parameters and dimensions
    lbm fl(Re,tau,D,nx,ny,odir);    

    // Call initialization functions to create simulation region and set up the initial condition
    fl.initialize(1.,0.,0.,flowtype);
    fl.solve(15000,100,bctype,forcetype,w,l);
    printf("Ran for 15000 iterations. Saved output every 150 frames. Finished!\n");

    // Run multiple passes for different widths w
    for (int w=10; w<=50; w += 10)
    {
        // Generate output directory for this choice of w
        char odir_w[50];
        sprintf(odir_w, "lbm_pbc_Re001.out/width_%d", w);

        // Create the simulation domain with the specified parameters and dimensions
        lbm fl(Re,tau,D,nx,ny,odir_w);

        // Call initialization functions to create simulation region and set up the initial condition
        fl.initialize(1.,0.,0.,flowtype);
        fl.solve(15000,100,bctype,forcetype,w,l);
        printf("Ran for 15000 iterations with w=%d. Saved output every 150 frames. Finished!\n", w);
    }
}