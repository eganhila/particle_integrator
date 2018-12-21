#ifndef GTEST
#include <iostream>
#include "sim_dat.h"
#include "sim_controller.h"
#include "integrator.h"
#include "acceleration.h"

int main(){


    SimDat sd(120);
    read_simulation_data(sd,
            "/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_10nT_Eint.h5");

    SimController sc(100, 10000, 0.01, sd, "pinter_output.h5");
    sc.set_particle_pop(16, 1, 1000000.0);
    sc.uniform_E = true;
    sc.run_cell(59,59, 39);
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
