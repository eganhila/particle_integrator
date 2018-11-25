#ifndef INTEGRATOR_H
#define INTEGRATOR_H

struct Particle{
        float state[6];
        float mass;     
        float charge;
};

struct Derivative {
        float d[6]={0};
};

void acceleration(const Particle & particle, double t, float * acc);
Derivative evaluate(const Particle &particle_init, double t, float dt, const Derivative & d);
void integrate(Particle & particle, double t, float dt);

#endif /* INTEGRATOR_H */
