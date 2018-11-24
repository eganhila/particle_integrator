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
        out[0] = Bx[i,j,k];
        out[1] = By[i,j,k];
        out[2] = Bz[i,j,k];
        out[3] = Ex[i,j,k];
        out[4] = Ey[i,j,k];
        out[5] = Ez[i,j,k];
    }


    int dim;
    float * Bx,* By, * Bz;
    float * Ex,* Ey, * Ez;
    float * x, * y, * z;
};

int read_simulation_data(SimDat& sd);
