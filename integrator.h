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
        float t0, t, t_final;
        float dt;
        Particle * particle;
        SimDat * sd = NULL;
        bool has_sd = false;
        void (*acc_func)(const Particle & particle, const SimDat &sd, float * acc) = NULL;


        void set_sd(SimDat & new_sd){ sd = &new_sd;
            has_sd = true;
            }
        
        void set_accel(void new_acc_func(const Particle & particle, const SimDat &sd, float * acc)){acc_func = new_acc_func;}


        Integrator(Particle & new_particle, float new_t0, float new_tfinal,  float new_dt) 
            : t0(new_t0), t(new_t0), t_final(new_tfinal),dt(new_dt), particle(&new_particle){}

        void evaluate_derivative(float t, float dt, const float * d_in,  float * d_out);
        bool integrate_step();
        int integrate();
        int evaluate_bcs();
};




#endif /* INTEGRATOR_H */
