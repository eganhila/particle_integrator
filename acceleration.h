#ifndef ACCEL_H
#define ACCEL_H

#include <iostream>
#include "integrator.h"
#include "sim_dat.h"
#include "interpolate.h"


bool zero_accel(const Particle & particle, const SimDat &sd, float * acc){
    for (int i=0; i<3; i++){
        acc[i] = 0;}
    return true;
};

bool constBz_accel(const Particle & particle, const SimDat &sd, float * acc){
        const float B[3] = {0,0,2};

        for (int axis=0; axis<3; axis++){
            acc[axis] =  particle.charge *( -1*particle.state[3+(axis+1)%3]*B[(axis+2)%3]+particle.state[3+(axis+2)%3]*B[(axis+1)%3])/particle.mass;
        }
    return true;
};

bool simDat_accel(const Particle & particle, const  SimDat &sd, float * acc){
    float pss[6]; 
    bool status;

    //Get interpolated data
    status = TrilinearInterpolate(particle.state, sd, pss);

    //Calc acceleration
    for (int i=0; i<3; i++){
        //factor for Z/mu and m -> km
        acc[i] = 0.09572*particle.charge * (
                -1*particle.state[3+(i+1)%3]*pss[(i+2)%3]+
                particle.state[3+(i+2)%3]*pss[(i+1)%3]+pss[3+i])/particle.mass;
    }
    return status;

};


#endif
