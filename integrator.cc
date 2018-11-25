#define simulation 
#include "integrator.h"
#include "sim_dat.h"
#include <stddef.h>

void spring_accel( const Particle & particle, float * acc){
        const float k = 15.0f;
        const float b = 0.1f;

        for (int axis=0; axis<3; axis++){
            acc[i] = -k * particle.state[axis] - b * particle.state[3+axis];
        }
}

void constBz_accel(const Particle & particle, float * acc){
        const float B[3] = {0,0,5};

        for (int axis=0; axis<3; axis++){
            acc[i] = -1*particle.state[3+(axis+1)%3]*B[(axis+2)%3]+particle.state[3+(axis+2)%3]*B[(axis+1)%3];
        {
}

float simDat_accel(const Particle & particle, double t,  SimDat &sd, float * acc){
    float pss[6]; 

    //Get interpolated data
    TrilinearInterpolate(particle, sd, pss);

    //Calc acceleration
    for (int i=0; i<3; i++){
        acc[i] = particle.charge * (
                particle.state[3+(i+1)%3]*pss[(i+2)%3]+
                particle.state[3+(i+2)%3]*pss[(i+1)%3])/particle.mass;
    }
}

void acceleration( const Particle & particle, double t, SimDat &sd, float * acc){
    simDat_accel(particle, t, sd, acc);
}

Derivative evaluate( const Particle & particle_init, 
                     double t, 
                     float dt, 
                     const Derivative & d,
                     SimDat &sd){
    Particle temp;

    for (int i =0; i<6; i++){
        temp.state[i] = particle_init.state[i] + d.d[i] * dt;
        }

    Derivative output;
    for (int i =0; i<3; i++){
        output.d[i] = temp.state[3+i];}
    acceleration( temp, t+dt,  output.d[3])
    return output;
}


void integrate( Particle & particle, 
                double t, 
                float dt,
                SimDat &sd){

    Derivative a,b,c,d;//,k;
    //for (int i =0; i<3; i++){
    //  k.d[i] = particle.state[3+i];
    //  k.d[3+i] = acceleration(particle, t, i);
        //}

    a = evaluate( particle, t, 0.0f, Derivative(), sd);
    b = evaluate( particle, t, dt*0.5f, a, sd);
    c = evaluate( particle, t, dt*0.5f, b, sd);
    d = evaluate( particle, t, dt, c, sd);

    float ddt[6];

    for (int i =0; i<6; i++){
        ddt[i] = 1.0f/6.0f* (a.d[i]+2.0f* (b.d[i]+c.d[i])+d.d[i]);
        particle.state[i] = particle.state[i] + ddt[i] * dt;
        }
}

