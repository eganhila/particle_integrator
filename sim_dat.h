
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

int read_simulation_data(SimDat& sd);
