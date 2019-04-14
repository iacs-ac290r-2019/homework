#include <cstdio>

#include "lbm.hh"

int main() {
    char odir[]="lbm_mbc_Re001.out"; // output directory

    /** This is a simple Poiseuille flow.

     *  Re = 0.01
     *  nu = 0.1 -> tau = 0.8
     *  umax = Re*nu/D = 0.01*0.1/60 = 0.001/60
     *  mu = rho*nu
     *  f = 8*mu*umax/H/H = 8*0.1*0.001/60/60/60

     *  Periodic boundary condition on west and east.
     *  No-slip on north and south.
     *  initial condition: rho_0 = 1.0, ux_0 = 0.0, uy_0 = 0.0
     *  Expected outcome: uy stays at 0.0, density stays at 1.0,
     *                    ux has parabolic behavior at vertical cut.
     *  forcing gives the velocity. 
     *  Run 10000 time (relative) and output every 100 frame. */

    // Set up simulation parameters.
    double Re=0.01;      // Reynolds number
    double tau=0.8;    // Relaxation constant
    double D=60;       // Reference length    
    // Set up grid properties.
    int nx=200;        // Channel length
    int ny=60;         // Channel width 
    int w=30;          // Narrowing width
    int l=50;          // Narrowing length
    int flowtype=1;    // 0: Standard, 1: Inlet-Outlet
    int bctype=1;      // 0: Periodic, 1: Open
    double forcetype=0.;
    
    // Create the simulation domain with the specified parameters and dimensions
    lbm fl(Re,tau,D,nx,ny,odir);

    // Calculate the input velocity
    // double u0=Re*fl.nu/D;
    // printf("Sanity check u0: %g\n",u0);

    // Call initialization functions to create simulation region and set up the initial condition
    fl.initialize(1.,0.,0.,flowtype);
    // double forcetype=0.0001;
    fl.solve(1000,100,bctype,forcetype,w,l);
    printf("Ran for 2000 iterations. Saved output every 20 frames. Finished!\n");

    // Run multiple passes for different widths w
    for (int w=10; w<=50; w += 10)
    {
        // Generate output directory for this choice of w
        char odir[50];
        sprintf(odir, "lbm_mbc_Re001.out/width_%d", w);

        // Create the simulation domain with the specified parameters and dimensions
        lbm fl(Re,tau,D,nx,ny,odir);

        // Call initialization functions to create simulation region and set up the initial condition
        fl.initialize(1.,0.,0.,flowtype);
        fl.solve(2000,100,bctype,forcetype,w,l);
        printf("Ran for 2000 iterations with w=%d. Saved output every 20 frames. Finished!\n", w);
    }
}