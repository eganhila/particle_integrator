#ifndef GTEST
#include <iostream>
#include "sim_dat.h"
#include "sim_controller.h"
#include "integrator.h"

int main(){


    SimDat sd(120);
    read_simulation_data(sd,
            "/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_50nT_Eint.h5");


    SimController sc(10, 100, 0.01, sd);
    sc.set_particle_pop(16, 1, 200);
    sc.run();


    return 0;
}
#endif
