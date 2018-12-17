#ifndef GTEST
#include <iostream>
#include "sim_dat.h"
#include "sim_controller.h"
#include "integrator.h"
#include "acceleration.h"

int main(){


    SimDat sd(120);
    read_simulation_data(sd,
            "/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_100nT_Eint.h5");


    SimController sc(8, 1000, 0.01, sd, "pinter_output.h5");
    sc.set_particle_pop(16, 1, 300.0);
    sc.run();

    /*
    Particle p;
    p.state[0] = sd.x[60];
    p.state[1] = sd.y[40];
    p.state[2] = sd.z[60];
    p.state[3] = 0;
    p.state[4] = 0;
    p.state[5] = 0;

    Integrator intg(p, 1,1000, 0.1);
    intg.set_sd(sd);
    intg.set_accel(simDat_accel);
    intg.verbose = true;
            
    //Particle integrated (work done here)
    int status;
    status = intg.integrate();
    //std::cout<<status<<std::endl;
    */
    return 0;
}
#endif
