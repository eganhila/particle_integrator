#define const_Bz
#include <iostream>
#include "hdf5.h"

#define FILE            "/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_0nT.h5"
#define DATASET     "magnetic_field_x"
#define DIM            120
using namespace std;

struct SimDat {

    SimDat(int new_dim){
        dim = new_dim;
        int dim3 = dim*dim*dim;

        x = new float[dim];
        y = new float[dim];
        z = new float[dim];
        Bx = new float[dim3];
        By = new float[dim3];
        Bz = new float[dim3];
        Ux = new float[dim3];
        Uy = new float[dim3];
        Uz = new float[dim3];

    }

    ~SimDat(){
        delete[] x;
        delete[] y;
        delete[] z;
        delete[] Bx;
        delete[] By;
        delete[] Bz;
        delete[] Ux;
        delete[] Uy;
        delete[] Uz;
    }

    int dim;
    float * Bx,* By, * Bz;
    float * Ux,* Uy, * Uz;
    float * x, * y, * z;
};


int read_simulation_data(SimDat& sd){
    hid_t       file, space, dset, dsetx, dsety, dsetz;          /* Handles */
    herr_t      status;
    hsize_t     dims[3] = {DIM, DIM, DIM};
    float         rdata[DIM][DIM][DIM];          /* Read buffer */
    int            i, j, k;

    file = H5Fopen (FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

    //x,y,z
    dsetx = H5Dopen (file, "x", H5P_DEFAULT);
    dsety = H5Dopen (file, "y", H5P_DEFAULT);
    dsetz = H5Dopen (file, "z", H5P_DEFAULT);

    status = H5Dread(dsetx, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.x);
    status = H5Dread(dsety, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.y);
    status = H5Dread(dsetz, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.z);

    status = H5Dclose (dsetx);
    status = H5Dclose (dsety);
    status = H5Dclose (dsetz);
    
    //mag field
    dsetx = H5Dopen (file, "magnetic_field_x", H5P_DEFAULT);
    dsety = H5Dopen (file, "magnetic_field_y", H5P_DEFAULT);
    dsetz = H5Dopen (file, "magnetic_field_z", H5P_DEFAULT);

    status = H5Dread(dsetx, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Bx);
    status = H5Dread(dsety, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.By);
    status = H5Dread(dsetz, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Bz);

    status = H5Dclose (dsetx);
    status = H5Dclose (dsety);
    status = H5Dclose (dsetz);
    
    //v field
    dsetx = H5Dopen (file, "electron_velocity_x", H5P_DEFAULT);
    dsety = H5Dopen (file, "electron_velocity_y", H5P_DEFAULT);
    dsetz = H5Dopen (file, "electron_velocity_z", H5P_DEFAULT);

    status = H5Dread(dsetx, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Ux);
    status = H5Dread(dsety, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Uy);
    status = H5Dread(dsetz, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Uz);

    status = H5Dclose (dsetx);
    status = H5Dclose (dsety);
    status = H5Dclose (dsetz);

    status = H5Fclose (file);
}


struct Particle{
        float state[6];
        float mass;     // velocity
};

struct Derivative {
        float d[6]={0};
};


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


int main(){
        Particle particle;
        SimDat sd(120);
        float t=0, dt=0.001, T=2;

        particle.state[0] = 0;
        particle.state[1] = 0;
        particle.state[2] = 0;
        particle.state[3] = 1;
        particle.state[4] = 1;
        particle.state[5] = 1;

        read_simulation_data(sd);
/*
        while (t < T){
                integrate(particle, t, dt);
                std::cout << particle.state[0]<<" " << particle.state[1]<<"\n";
                t += dt;
        }
*/        


        return 0;
}
