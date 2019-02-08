#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "sim_dat.h"
#include <iostream>
struct Particle{
        float state[6];
        float mass=1;     
        float charge=1;

        void print(){
            std::cout<<"Particle: ["<<state[0]<<", "<<state[1];
            std::cout<<", "<<state[2]<<", "<<state[3];
            std::cout<<", "<<state[4]<<", "<<state[5]<<"] ";
            std::cout<<", q: "<<charge<<", mu: "<<mass<<std::endl;
        }
        void print_state(){
            std::cout<<state[0]<<", "<<state[1];
            std::cout<<", "<<state[2]<<", "<<state[3];
            std::cout<<", "<<state[4]<<", "<<state[5]<<std::endl;
        }
};

class Integrator {
    public:
        float t0, t, t_final;
        float dt;
        bool verbose=false;
        Particle * particle;
        SimDat * sd = NULL;
        bool has_sd = false;
        float init_state[6];
        bool within_psphere = false;
        bool within_tail = false;
        bool (*acc_func)(const Particle & particle, const SimDat &sd, float * acc) = NULL;


        void set_sd(SimDat & new_sd){ sd = &new_sd;
            has_sd = true;
            }
        
        void set_accel(bool new_acc_func(const Particle & particle, const SimDat &sd, float * acc)){acc_func = new_acc_func;}


        Integrator(Particle & new_particle, float new_t0, float new_tfinal,  float new_dt) 
            : t0(new_t0), t(new_t0), t_final(new_tfinal),dt(new_dt), particle(&new_particle){}

        bool evaluate_derivative(float t, float dt, const float * d_in,  float * d_out);
        bool integrate_step();
        int integrate();
        int evaluate_bcs();
        int evaluate_final_status();
        bool evaluate_tail();
        bool evaluate_psphere();
};




#endif /* INTEGRATOR_H */
