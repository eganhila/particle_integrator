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
    /*
    for (int i =0; i<3; i++){
        ddt[i] = 1.0f/6.0f* (d1[i]+2.0f* (d2[i]+d3[i])+d4[i]);
        std::cout << i<<": "<< particle->state[i]<<", "<< ddt[i] <<", "<< particle->state[i+3]<<", "<<ddt[i+3]<<std::endl; 
    }*/

    for (int i =0; i<6; i++){
        ddt[i] = 1.0f/6.0f* (d1[i]+2.0f* (d2[i]+d3[i])+d4[i]);
        particle->state[i] = particle->state[i] + ddt[i];
        }
    t = t+dt;


    //std::cout << particle->state[0]<<" " <<particle->state[1]<<" " <<particle->state[2]<<" " <<particle->state[3]<<" " << particle->state[4]<< " "<<particle->state[5] <<"\n";
    return success;
}

int Integrator::evaluate_bcs(){
    float planet_r2 = 3390.0*3390.0, particle_r2;

    for (int axis=0; axis<3; axis++){
        if ((particle->state[axis]+2*particle->state[3+axis]*dt < sd->bbox[axis]) || 
            (particle->state[axis]+2*particle->state[3+axis]*dt > sd->bbox[axis+3])){
            return 1;
            }
    }
    particle_r2 = particle->state[0]*particle->state[0]+
                  particle->state[1]*particle->state[1]+ 
                  particle->state[2]*particle->state[2]; 

    if (particle_r2 < planet_r2){return 2;}

    return 0;

}


int Integrator::integrate(){
    int pstatus = 0;

    while((t<t_final)){
        integrate_step();
        
        //check if still in bounds
        if (has_sd) {
           pstatus = evaluate_bcs(); 
           if (pstatus != 0){ break;}
        }
    }
    return pstatus;
}
