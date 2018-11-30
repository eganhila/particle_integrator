#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "sim_dat.h"
struct Particle{
        float state[6];
        float mass;     
        float charge;
};

void acceleration(const Particle & particle, double t, const SimDat &sd, float * acc);
void evaluate_derivative(const Particle &particle_init, double t, float dt, const float * d_in, const SimDat &sd, float * d_out);
bool integrate(Particle & particle, double t, float dt, SimDat &sd);

#endif /* INTEGRATOR_H */
