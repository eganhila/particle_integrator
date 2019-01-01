#define simulation 
#include "integrator.h"
#include "sim_dat.h"
#include "interpolate.h"
#include <stddef.h>
#include <iostream>
#include "math.h"



bool Integrator::evaluate_derivative( 
                     float t, 
                     float dt, 
                     const float * d_in,
                     float * d_out){
    float * y_n = particle->state;
    bool success;
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
    success = acc_func(temp, *sd,  &d_out[3]);
    if (!success){return false;}

    for (int i=0;i<3;i++){
        d_out[3+i] = d_out[3+i]*dt;
    }

    return true;
}


bool Integrator::integrate_step(){

    bool success = true;
    float  d0[6] = {0,0,0,0,0,0},d1[6],d2[6],d3[6], d4[6], temp[6];

    //This needs to be float checked b/c edited eval deriv func
    success = evaluate_derivative( t, dt, d0, d1);
    if (!success) {return false;}
    
    for (int i=0; i<6; i++){temp[i] = d1[i]/2.0;}
    success = evaluate_derivative( t+dt*0.5f, dt, temp,  d2);
    if (!success) {return false;}

    for (int i=0; i<6; i++){temp[i] = d2[i]/2.0;}
    success = evaluate_derivative( t+dt*0.5f, dt, temp, d3);
    if (!success) {return false;}

    success = evaluate_derivative( t+dt, dt, d3,  d4);
    if (!success) {return false;}

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

    if (verbose){particle->print_state();}

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

bool Integrator::evaluate_psphere(){
    if (particle->state[0]>0.5*3390+init_state[0]){
        return true;}
    else {return false;}

}
bool Integrator::evaluate_tail(){
    if (particle->state[0]<-1.1*3390){
        return true;}
    else {return false;}

}

int Integrator::evaluate_final_status(){
    float vx,vy,vz,vfinal,dr_x,dr_y,dr_z,dr;
    int pstatus;
    pstatus = evaluate_bcs();

    //escape
    if (pstatus == 1){
        if (within_psphere){
            return 3; //plasmasphere escape
        }
        else{
            vx = particle->state[3]; 
            vy = particle->state[4]; 
            vz = particle->state[5]; 
            vfinal = pow(vx*vx+vy*vy+vz*vz, 0.5);

            if (vfinal > 0.95*350){
                return 6; //plume escape
            }
            else {
                return 5; //Direct escape
            }
        }
    }
    else{ //bound

        dr_x = particle->state[0]-init_state[0];
        dr_y = particle->state[1]-init_state[1];
        dr_z = particle->state[2]-init_state[2];

        dr = pow(dr_x*dr_x+dr_y*dr_y+dr_z*dr_z,0.5);

        if (dr < 0.5*3390){
            return 1; //closely bound
        }
        else{

            if (within_tail){ return 2;} //bound tail
            else {return 3;} //bound pshere

        }

    }
}


int Integrator::integrate(){
    int pstatus = 0;
    bool success;


    while((t<t_final)){
        success = integrate_step();
        if (!success){
            pstatus = -1;
            break;
        }
        
        //check if still in bounds
        if (has_sd) {
           pstatus = evaluate_bcs(); 
           if (! within_psphere){
               within_psphere = evaluate_psphere();}
           if (! within_tail){
               within_tail = evaluate_tail();
           }
           if (pstatus != 0){ break;}
        }
    }

    pstatus = evaluate_final_status();
    return pstatus;
}
