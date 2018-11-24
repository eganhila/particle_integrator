#define const_Bz
#include "integrator.h"

float acceleration( const Particle & particle, double t, int axis){

    #ifdef spring
        const float k = 15.0f;
        const float b = 0.1f;

        return -k * particle.state[axis] - b * particle.state[3+axis];
    #endif

     #ifdef const_Bz
        const float B[3] = {0,0,5};
        return -1*particle.state[3+(axis+1)%3]*B[(axis+2)%3]+particle.state[3+(axis+2)%3]*B[(axis+1)%3];
     #endif

}

Derivative evaluate( const Particle & particle_init, 
                     double t, 
                     float dt, 
                     const Derivative & d ){
    Particle temp;

    for (int i =0; i<6; i++){
        temp.state[i] = particle_init.state[i] + d.d[i] * dt;
        }

    Derivative output;
    for (int i =0; i<3; i++){
        output.d[i] = temp.state[3+i];
        output.d[3+i] = acceleration( temp, t+dt, i);
        }
    return output;
}


void integrate( Particle & particle, 
                double t, 
                float dt ){

    Derivative a,b,c,d;//,k;
    //for (int i =0; i<3; i++){
    //  k.d[i] = particle.state[3+i];
    //  k.d[3+i] = acceleration(particle, t, i);
        //}

    a = evaluate( particle, t, 0.0f, Derivative());
    b = evaluate( particle, t, dt*0.5f, a );
    c = evaluate( particle, t, dt*0.5f, b );
    d = evaluate( particle, t, dt, c );

    float ddt[6];

    for (int i =0; i<6; i++){
        ddt[i] = 1.0f/6.0f* (a.d[i]+2.0f* (b.d[i]+c.d[i])+d.d[i]);
        particle.state[i] = particle.state[i] + ddt[i] * dt;
        }
}

