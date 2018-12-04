#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "sim_dat.h"
struct Particle{
        float state[6];
        float mass=1;     
        float charge=1;
};

class Integrator {
    public:
        double t0, t;
        float dt;
        Particle * particle;
        SimDat * sd = NULL;
        void (*acc_func)(const Particle & particle, const SimDat &sd, float * acc) = NULL;


        void set_sd(SimDat & new_sd){ sd = &new_sd;
            }
        
        void set_accel(void new_acc_func(const Particle & particle, const SimDat &sd, float * acc)){acc_func = new_acc_func;}


        Integrator(Particle & new_particle, double new_t0, float new_dt) 
            : t0(new_t0), t(new_t0), dt(new_dt), particle(&new_particle){}

        void evaluate_derivative(double t, float dt, const float * d_in,  float * d_out);
        bool integrate_step();
        bool integrate();
};




#endif /* INTEGRATOR_H */
