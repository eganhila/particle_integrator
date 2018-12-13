#include "sim_dat.h"
#include "integrator.h"
#include "sim_controller.h"
#include "acceleration.h"
#include <random>
#include <iostream>
#include "hdf5.h"
#include <math.h>


float inv_maxwell_cdf(float v, float a){
    return erf(v/(pow(2,0.5)*a))- pow(2/M_PI, 0.5)*(v*exp(-1*v*v/(2*a*a)))/a;

}


Particle SimController :: draw_particle(int cell_idx){

    float x, y, z, x0,y0,z0,vx,vy,vz,u,v;
    float noon_f=1.0, night_f=0.1;

    // Get base position
    x0 = sd->x[cell_idx];
    y0 = sd->y[cell_idx];
    z0 = sd->z[cell_idx];

    //setup randoming
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> uniform_dist(0, 1);

    //Spatially distributed uniformly in cell
    x = x0 + (uniform_dist(e2)-0.5)*sd->dx;
    y = y0 + (uniform_dist(e2)-0.5)*sd->dx;
    z = z0 + (uniform_dist(e2)-0.5)*sd->dx;


    //draw 3 random components and normalize to
    //  get uniformly distributed normalized
    //  velocity vector
    vx = uniform_dist(e2);
    vy = uniform_dist(e2);
    vz = uniform_dist(e2);

    v = pow(vx*vx+vy*vy+vz*vz,0.5);
    vx = vx/v;
    vy = vy/v;
    vz = vz/v;

    //Now want to draw speed from maxwellian
    //  distribution using inverse random sampling
    //  and scale velocity vector to maxwell speed
    //  90.82e-3 factor from sqrt(kT/mu mH) into km/s

    u = uniform_dist(e2);
    v = inv_maxwell_cdf(u, 90.82e-3*pow(pop_T/pop_m,0.5));

    vx = vx*v;
    vy = vy*v;
    vz = vx*v;


    Particle p;
    p.state[0] = x;
    p.state[1] = y;
    p.state[2] = z;
    p.state[3] = vx;
    p.state[4] = vy;
    p.state[5] = vz;
    p.mass = pop_m; 
    p.charge = pop_q;
    return p;
}

void SimController :: run(){
    int cell_idx;
    int p_idx;
    Particle * cell_particles;
    int * all_status;
    float * positions, * velocities;
    hid_t       file_id;   /* file identifier */
    herr_t      status;

    std::cout<<"setting up"<<std::endl;
    /* Create a new file using default properties. */
    setup_datawriter();

    // Setup data structures that we will
    // keep re-using
    cell_particles = new Particle[N_particles]; 
    all_status = new int[N_particles];
    positions = new float[N_particles*3];
    velocities = new float[N_particles*3];


    // Iterate over all cells
    for (cell_idx=0; cell_idx<sd->dim3; cell_idx++){

        //Evaluate if we want to check this cell
        if (!eval_cell(cell_idx)){continue;}
        

        // Each cell needs N_particles, this is what we're going
        // to parallelize over
        for (p_idx=0; p_idx< N_particles; p_idx++){

            Particle p;
            //draw particle from distribution
            p = draw_particle(cell_idx);
            cell_particles[p_idx]=p;

            //store init particle state
            for (int i=0;i<3;i++){
                positions[p_idx+i*N_particles] = p.state[i];
                velocities[p_idx+i*N_particles] = p.state[i+3];
            }


            //Setup integrator
            Integrator intg(p, 1,100, 0.1);
            intg.set_sd(*sd);
            intg.set_accel(simDat_accel);
            
            //Particle integrated (work done here)
            all_status[p_idx] = intg.integrate();
            
        }

        // Write out cell
        write_cell_data(cell_idx, positions, velocities, all_status);
    } 

    //Delete data structures
    delete[] cell_particles;
    delete[] all_status;
    delete[] positions; 
    delete[] velocities;



}

void SimController :: setup_datawriter(){
    hid_t       file_id, dataset_id, dataspace_id, group;  /* identifiers */
    hsize_t  dims2D[2], dims3D[3];
    herr_t      status;
    char sidx[20], dset_name[30];

    dims2D[0] = sd->dim3;
    dims2D[1] = N_particles;

    dims3D[0] = 3;
    dims3D[1] = sd->dim3;
    dims3D[2] = N_particles;

    file_id = H5Fcreate(outname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dataspace_id = H5Screate_simple(2, dims2D, NULL);
    dataset_id = H5Dcreate2(file_id, "status", H5T_STD_I32BE,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    dataspace_id = H5Screate_simple(3, dims3D, NULL);
    dataset_id = H5Dcreate2(file_id, "position", H5T_NATIVE_FLOAT,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(3, dims3D, NULL);
    dataset_id = H5Dcreate2(file_id, "velocity", H5T_NATIVE_FLOAT,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    status = H5Fclose(file_id);
}

void SimController :: write_cell_data(int cell_idx, float * positions, float * velocities, int * all_status){
    hid_t       file_id, dataset_id, dataspace_id, memspace_id;  /* identifiers */
    hsize_t  dims2D[2], dims3D[3], offset2D[2], offset3D[3];
    hsize_t stride[3] = {1,1,1}, block[3]={1,1,1};
    herr_t      status;
    char sidx[20], dset_name[30];

    dims2D[0] = 1;
    dims2D[1] = N_particles;

    dims3D[0] = 3;
    dims3D[1] = 1;
    dims3D[2] = N_particles;

    offset2D[0] = cell_idx;
    offset2D[1] = 0;
    offset3D[0] = 0;
    offset3D[1] = cell_idx;
    offset3D[2] = 0;



    // open existing file
    file_id = H5Fopen(outname, H5F_ACC_RDWR, H5P_DEFAULT);

    //create/close dataspace, dataset, write data, for status
    dataset_id = H5Dopen2(file_id, "status", H5P_DEFAULT);
    memspace_id = H5Screate_simple(2,dims2D, NULL); 
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D,
                                 stride, dims2D, block);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
                      all_status);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    //create/close dataspace, dataset, write data, for init location
    dataset_id = H5Dopen2(file_id, "position", H5P_DEFAULT);
    memspace_id = H5Screate_simple(3,dims3D, NULL); 
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D,
                                 stride, dims3D, block);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
                      positions);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    //create/close dataspace, dataset, write data, for init velocity
    dataset_id = H5Dopen2(file_id, "velocity", H5P_DEFAULT);
    memspace_id = H5Screate_simple(3,dims3D, NULL); 
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D,
                                 stride, dims3D, block);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
                      velocities);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    //close group and file
    status = H5Fclose(file_id);
}

bool SimController :: eval_cell(int cell_idx){
    bool evaluate = true;
    float x,y,z,r;
    x = sd->x[cell_idx];
    y = sd->y[cell_idx];
    z = sd->z[cell_idx];
    r = pow(x*x+y*y+z*z,0.5)/3390.0;

    //want < 2.5 RM, > 1 RM
//    if (r > 1.5){evaluate=false;}
//    if (r < 1){evaluate=false;}

    evaluate=false;
    if (cell_idx==23979){evaluate=true;}
    if (cell_idx==223749){evaluate=true;}
    
    return evaluate;
}

void SimController :: set_particle_pop(float mass, float charge, float temperature){
    pop_m = mass;
    pop_q = charge;
    pop_T = temperature;

    pop_set=true;
}
