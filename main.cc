#include <iostream>
#include "sim_dat.h"
#include "integrator.h"

#ifndef GTEST
int main(){
        std::cout <<"Running Main"<<std::endl; 
        Particle particle;
        SimDat sd(120);
        float t=0, dt=0.001, T=2;
        bool success;

        particle.state[0] = 1.5;
        particle.state[1] = 0;
        particle.state[2] = 0;
        particle.state[3] = 1;
        particle.state[4] = 1;
        particle.state[5] = 1;

        read_simulation_data(sd, "/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_0nT.h5");

        Integrator integrator(particle, t0, dt);
        integrator.set_accel(constBz_accel);


        while (t < T){
                integrator.integrate(particle, t, dt, sd);
                std::cout << particle.state[0]<<" " << particle.state[1]<<"\n";
                t += dt;
                if (success == false){ break;}
        }
        


        return 0;
}
#endif
