#define simulation 
#include "integrator.h"
#include "sim_dat.h"
#include "interpolate.h"
#include <stddef.h>
#include <iostream>



void Integrator::evaluate_derivative( 
                     float t, 
                     float dt, 
                     const float * d_in,
                     float * d_out){
    float * y_n = particle->state;
    Particle temp;

    //make fake particle at propper time
    for (int i =0; i<6; i++){
        temp.state[i] = particle->state[i] + d_in[i];
    }
    temp.mass = particle->mass; 
    temp.charge = particle->charge;

    // spatial output
    for (int i =0; i<3; i++){
        d_out[i] = dt*temp.state[3+i]; 
    }


    //Use fake particle to calc velocity space deriv
    acc_func(temp, *sd,  &d_out[3]);

    for (int i=0;i<3;i++){
        d_out[3+i] = d_out[3+i]*dt;
    }
}


bool Integrator::integrate_step(){

    bool success = true;
    float  d0[6] = {0,0,0,0,0,0},d1[6],d2[6],d3[6], d4[6], temp[6];

    //This needs to be float checked b/c edited eval deriv func
    evaluate_derivative( t, dt, d0, d1);
    
    for (int i=0; i<6; i++){temp[i] = d1[i]/2.0;}
    evaluate_derivative( t+dt*0.5f, dt, temp,  d2);

    for (int i=0; i<6; i++){temp[i] = d2[i]/2.0;}
    evaluate_derivative( t+dt*0.5f, dt, temp, d3);

    evaluate_derivative( t+dt, dt, d3,  d4);

    float ddt[6];

    for (int i =0; i<6; i++){
        ddt[i] = 1.0f/6.0f* (d1[i]+2.0f* (d2[i]+d3[i])+d4[i]);
        particle->state[i] = particle->state[i] + ddt[i];
        }
    t = t+dt;

    //std::cout << particle->state[0]<<" " <<particle->state[1]<<" " <<particle->state[2]<<" " <<particle->state[3]<<" " << particle->state[4]<< " "<<particle->state[5] <<"\n";
    return success;
}


bool Integrator::integrate(){
    bool success = true;
    bool within_box = true;

    while((t<t_final)&& within_box){
        integrate_step();
        
        //check if still in bounds
        if (has_sd) {
            
            for (int axis=0; axis<3; axis++){
                if ((particle->state[axis]+2*particle->state[3+axis]*dt < sd->bbox[axis]) || 
                    (particle->state[axis]+2*particle->state[3+axis]*dt > sd->bbox[axis+3])){
                    within_box = false;
                    }
            }
        }
    }

    return success;
}
