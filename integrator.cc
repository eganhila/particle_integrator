#define simulation 
#include "integrator.h"
#include "sim_dat.h"
#include "interpolate.h"
#include <stddef.h>

void spring_accel( const Particle & particle, float * acc){
        const float k = 15.0f;
        const float b = 0.1f;

        for (int axis=0; axis<3; axis++){
            acc[axis] = -k * particle.state[axis] - b * particle.state[3+axis];
        }
}

void constBz_accel(const Particle & particle, float * acc){
        const float B[3] = {0,0,5};

        for (int axis=0; axis<3; axis++){
            acc[axis] = -1*particle.state[3+(axis+1)%3]*B[(axis+2)%3]+particle.state[3+(axis+2)%3]*B[(axis+1)%3];
        }
}

void simDat_accel(const Particle & particle, double t, const  SimDat &sd, float * acc){
    float pss[6]; 

    //Get interpolated data
    TrilinearInterpolate(particle.state, sd, pss);

    //Calc acceleration
    for (int i=0; i<3; i++){
        acc[i] = particle.charge * (
                particle.state[3+(i+1)%3]*pss[(i+2)%3]+
                particle.state[3+(i+2)%3]*pss[(i+1)%3])/particle.mass;
    }
}

void acceleration( const Particle & particle, double t, const SimDat &sd, float * acc){
    simDat_accel(particle, t, sd, acc);
}

void evaluate_derivative( const Particle & particle_init, 
                     double t, 
                     float dt, 
                     const float * d_in,
                     const SimDat &sd,
                     float * d_out){
    Particle temp;

    for (int i =0; i<6; i++){
        temp.state[i] = particle_init.state[i] + d_in[i] * dt;}

    for (int i =0; i<3; i++){
        d_out[i] = temp.state[3+i];}
    acceleration( temp, t+dt, sd,  &d_out[3]);
}


bool integrate( Particle & particle, 
                double t, 
                float dt,
                SimDat &sd){

    bool success = true;
    float  d0[6] = {0,0,0,0,0,0},d1[6],d2[6],d3[6], d4[6];//,k;
    //for (int i =0; i<3; i++){
    //  k.d[i] = particle.state[3+i];
    //  k.d[3+i] = acceleration(particle, t, i);
        //}

    evaluate_derivative( particle, t, 0.0f, d0, sd, d1);
    evaluate_derivative( particle, t, dt*0.5f, d1, sd, d2);
    evaluate_derivative( particle, t, dt*0.5f, d2, sd, d3);
    evaluate_derivative( particle, t, dt, d3, sd, d4);

    float ddt[6];

    for (int i =0; i<6; i++){
        ddt[i] = 1.0f/6.0f* (d1[i]+2.0f* (d2[i]+d3[i])+d4[i]);
        particle.state[i] = particle.state[i] + ddt[i] * dt;
        }

    return success;
}

