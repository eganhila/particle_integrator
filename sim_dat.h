#ifndef SIMDAT_H
#define SIMDAT_H

#include <string>

//Global Simulation Data
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
        Ex = new float[dim3];
        Ey = new float[dim3];
        Ez = new float[dim3];

    }

    ~SimDat(){
        delete[] x;
        delete[] y;
        delete[] z;
        delete[] Bx;
        delete[] By;
        delete[] Bz;
        delete[] Ex;
        delete[] Ey;
        delete[] Ez;
    }

    void GetSimState (const int i, const int j, const int k, float *out) const{
        int idx = i*dim*dim+j*dim+k;
        out[0] = Bx[idx]; 
        out[1] = By[idx];
        out[2] = Bz[idx];
        out[3] = Ex[idx];
        out[4] = Ey[idx];
        out[5] = Ez[idx];
    }


    int dim;
    float * Bx,* By, * Bz;
    float * Ex,* Ey, * Ez;
    float * x, * y, * z;
    float bbox[6];

    void set_bounds();
};

int read_simulation_data(SimDat& sd, const char* file_name);

#endif /* SIMDAT_H */
