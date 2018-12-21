#ifndef GTEST
#include <iostream>
#include "sim_dat.h"
#include "sim_controller.h"
#include "integrator.h"
#include "acceleration.h"
#include "math.h"
#include "interpolate.h"

int main(){


    SimDat sd(120);
    read_simulation_data(sd,
            "/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_10nT_Eint.h5");

    SimController sc(100, 10000, 0.01, sd, "pinter_output.h5");

    sc.set_particle_pop(16, 1, 1000000.0);
    sc.uniform_E = true;
    sc.setup_particlewriter();

//    sc.run_cell(59,59,39);

    float pos[3];
    int  idx[3];
    for (float theta=0; theta<=2*M_PI; theta+=0.1){
        pos[0] = cos(theta)*3690.0;
        pos[1] = 0;
        pos[2] = sin(theta)*3690.0;

        getCellIdx(sd, pos, idx);
        sc.run_cell(idx[0], idx[1], idx[2]);
    }
    
    //sc.run_sim();
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
}
#endif
