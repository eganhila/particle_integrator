#ifndef GTEST
#include <iostream>
#include "sim_dat.h"
#include "integrator.h"
#include "acceleration.h"

void setupConstSD(SimDat & sd){
    int d = sd.dim;
    int d3 = d*d*d;
        for (int i=0; i<(d3); i++){
            sd.Bx[i] = 0;
            sd.By[i] = 0;
            sd.Bz[i] = 2;
            sd.Ex[i] = 0;
            sd.Ey[i] = 0;
            sd.Ez[i] = 0;
        }
        for (int i=0; i<(d); i++){
            sd.x[i] = -4+8.0*i/d;
            sd.y[i] = -4+8.0*i/d;
            sd.x[i] = -4+8.0*i/d;
        }
        sd.set_bounds();
}
int main(){

    Particle p;
    p.state[0] = 2.5*3390;
    p.state[1] = 0;
    p.state[2] = 0;
    p.state[3] = -1;
    p.state[4] = 0;
    p.state[5] = 0;

    Integrator intg(p, 1,5, 0.01);

    SimDat sd(120);
    //setupConstSD(sd);
    read_simulation_data(sd, "/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_0nT.h5");
    intg.set_sd(sd);

    intg.set_accel(simDat_accel);
    intg.integrate();



    return 0;
}
#endif
