#ifndef GTEST
#include <iostream>
#include "sim_dat.h"
#include "sim_controller.h"
#include "integrator.h"
#include "acceleration.h"
#include "math.h"
#include "interpolate.h"

int main(){


    // Set filename and dimensionality of dataset to integrate through
    SimDat sd(120);
    read_simulation_data(sd,
            "/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_50nT_Eint.h5");

    /* Option 1: Integrate a single particle */

    /*
    Particle p;
    p.state[0] = 0*3390;
    p.state[1] = 1.3*3390;
    p.state[2] = -0*3390.0;
    p.state[3] = 0.809;
    p.state[4] = 0.199;
    p.state[5] = -0.77;
    p.mass = 16;

    Integrator intg(p, 0,1000, 0.01);
    intg.set_sd(sd);
    intg.set_accel(simDat_accel);
    intg.verbose = true;
            
    //Particle integrated (work done here)
    int status;
    status = intg.integrate();
    //std::cout<<status<<std::endl;
    return 0;
    */
    
    /* Option 2: Integrate a particle population */

    /*
    SimController sc(100, 100000, 0.1, sd, "pinter_output.h5");

    sc.set_particle_pop(16, 1, 10000.0);
    sc.uniform_E = false;
    sc.setup_particlewriter();
    */


    /* Can either run through a radius or a given cell(s) */
//    sc.run_radius(3390+300, 1e6);
//    sc.run_cell(59,59,39);

    /*
    float phi = M_PI/2; 
    float theta = 300;
    float r = 3839;
    float pos[3];
    int  idx[3];
    int i = 0;
    //for (float r=3390; r<=6780; r+=200){
        pos[0] = r*cos(theta)*sin(phi);
        pos[1] = r*sin(theta)*sin(phi); 
        pos[2] = r*cos(phi); 

        getCellIdx(sd, pos, idx);
        sc.run_cell(idx[0], idx[1], idx[2], i);
        i+=1;
   // }
*/


    // Or run with user defined boundaries in sim controller 
    //  sc.run_sim();
}
#endif
