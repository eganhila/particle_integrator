
#include "sim_dat.h"
#include "integrator.h"

int main(){
        Particle particle;
        SimDat sd(120);
        float t=0, dt=0.001, T=2;

        particle.state[0] = 0;
        particle.state[1] = 0;
        particle.state[2] = 0;
        particle.state[3] = 1;
        particle.state[4] = 1;
        particle.state[5] = 1;

        read_simulation_data(sd);
/*
        while (t < T){
                integrate(particle, t, dt);
                std::cout << particle.state[0]<<" " << particle.state[1]<<"\n";
                t += dt;
        }
*/        


        return 0;
}
