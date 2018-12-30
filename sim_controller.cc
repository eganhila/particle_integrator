#include "sim_dat.h"
#include "interpolate.h"
#include "integrator.h"
#include "sim_controller.h"
#include "acceleration.h"
#include <random>
#include <iostream>
#include "hdf5.h"
#include <math.h>


void reverse_cellidx(int cellidx, int dim, int & i, int &j,  int &k){
    i = cellidx/(dim*dim); 
    j = cellidx%(dim*dim)/dim;
    k = cellidx%(dim);
}


float maxwell_pdf(float v, float a){
    return pow(2/M_PI,0.5)*v*v*exp(-1*v*v/(2*a*a))/pow(a,3);
}


Particle SimController :: draw_radial_particle(float r){
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> uniform_dist(0, 1);

    Particle p = draw_particle(0);

    float xhat,yhat,zhat, mag;

    xhat = uniform_dist(e2)-0.5;
    yhat = uniform_dist(e2)-0.5;
    zhat = uniform_dist(e2)-0.5;

    mag = pow(xhat*xhat+yhat*yhat+zhat*zhat,0.5);

    xhat = xhat/mag;
    yhat = yhat/mag;
    zhat = zhat/mag;

    p.state[0] = xhat*r;
    p.state[1] = yhat*r;
    p.state[2] = zhat*r;

    return p;
}

Particle SimController :: draw_particle(int cell_idx){

    float x, y, z, x0,y0,z0,vx,vy,vz,pv,v, pu;
    float noon_f=1.0, night_f=0.1;
    int i,j,k;

    reverse_cellidx(cell_idx, sd->dim, i,j,k);

    // Get base position
    x0 = sd->x[i];
    y0 = sd->y[j];
    z0 = sd->z[k];

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


    float a = 90.82e-3*pow(pop_T/pop_m,0.5);
    if (uniform_E){
        v = uniform_dist(e2)*5*a;
    }
    else{
        //Changing to accept/reject method because before I 
        //  was using CDF not inverse CDF. 
        while (true){
            v = uniform_dist(e2)*5*a;
            pv = maxwell_pdf(v, a); 
            pu = uniform_dist(e2);
            if (pu<pv){break;}
        }
    }

    vx = vx*v;
    vy = vy*v;
    vz = vz*v;

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

void SimController :: run_radius(float r, int N_p){
    for (int p_idx=0; p_idx< N_p; p_idx++){

        Particle p;
        //draw particle from distribution
        p = draw_radial_particle(r);

        Integrator intg(p, 0.0, t_final, dt);
        intg.set_sd(*sd);
        intg.set_accel(simDat_accel);
        
        //Particle integrated (work done here)
        int Nsteps = 0;
        int p_status = 0;
        float out_data[6][MAX_STEPS] = {0}; 
        while ((p_status == 0)&&(Nsteps<MAX_STEPS)){
            intg.integrate_step();
            p_status = intg.evaluate_bcs();

            for (int i=0; i<6; i++){
                out_data[i][Nsteps] = p.state[i]; 
            }
            Nsteps +=1;
        }

        write_particle(0, p_idx, Nsteps, out_data);
    }
}

void SimController :: run_cell(int cell_i, int cell_j, int cell_k, int label){

    int cell_idx = cell_i*sd->dim*sd->dim+cell_j*sd->dim+cell_k;

    //setup_particlewriter();
    for (int p_idx=0; p_idx< N_particles; p_idx++){

        Particle p;
        //draw particle from distribution
        p = draw_particle(cell_idx);

        //Setup integrator
        Integrator intg(p, 0.0, t_final, dt);
        intg.set_sd(*sd);
        intg.set_accel(simDat_accel);
        
        //Particle integrated (work done here)
        int Nsteps = 0;
        int p_status = 0;
        float out_data[6][MAX_STEPS] = {0}; 
        while ((p_status == 0)&&(Nsteps<MAX_STEPS)){
            intg.integrate_step();
            p_status = intg.evaluate_bcs();

            for (int i=0; i<6; i++){
                out_data[i][Nsteps] = p.state[i]; 
            }
            Nsteps +=1;
        }

        if (label == -1){
            write_particle(cell_idx, p_idx, Nsteps, out_data);
        }else{
            write_particle(label, p_idx, Nsteps, out_data);
        }
    }

}

void SimController::setup_particlewriter(){
    hid_t file_id;
    herr_t      status;

    file_id = H5Fcreate(outname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Fclose(file_id);

}
void SimController::write_particle(int cell_idx, int p_idx, int N, float  pdata[6][MAX_STEPS]){
    herr_t      status;
    hid_t file_id, dataspace_id, dataset_id;

    // open existing file
    file_id = H5Fopen(outname, H5F_ACC_RDWR, H5P_DEFAULT);

    hsize_t dims[2];
    dims[0] = 6;
    dims[1] = MAX_STEPS;

    char buffer [20];
    sprintf (buffer, "%08d_%08d",cell_idx, p_idx);

    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate2(file_id, buffer, H5T_NATIVE_FLOAT,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, pdata);

    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    //close group and file
    status = H5Fclose(file_id);
}


void SimController :: run_sim(){
    int cell_idx;
    int p_idx;
    Particle * cell_particles;
    int * all_status;
    float * positions_start, * velocities_start, *positions_end, *velocities_end;
    hid_t       file_id;   /* file identifier */
    herr_t      status;

    /* Create a new file using default properties. */
    setup_datawriter();

    // Setup data structures that we will
    // keep re-using
    cell_particles = new Particle[N_particles]; 
    all_status = new int[N_particles];
    positions_start = new float[N_particles*3];
    velocities_start = new float[N_particles*3];
    positions_end = new float[N_particles*3];
    velocities_end = new float[N_particles*3];


    // Iterate over all cells
    for (cell_idx=0; cell_idx<sd->dim3; cell_idx++){

        //Evaluate if we want to check this cell
        if (!eval_cell(cell_idx)){continue;}
        

        // Each cell needs N_particles, this is what we're going
        // to parallelize over
        
        #pragma omp parallel for schedule(dynamic)
        for (p_idx=0; p_idx< N_particles; p_idx++){

            Particle p;
            //draw particle from distribution
            p = draw_particle(cell_idx);
            cell_particles[p_idx]=p;

            //store init particle state
            for (int i=0;i<3;i++){
                positions_start[p_idx+i*N_particles] = p.state[i];
                velocities_start[p_idx+i*N_particles] = p.state[i+3];
            }


            //Setup integrator
            Integrator intg(p, 0.0, t_final, dt);
            intg.set_sd(*sd);
            intg.set_accel(simDat_accel);
            
            //Particle integrated (work done here)
            all_status[p_idx] = intg.integrate();

            //store final particle state
            for (int i=0;i<3;i++){
                positions_end[p_idx+i*N_particles] = p.state[i];
                velocities_end[p_idx+i*N_particles] = p.state[i+3];
            }

        }
            

        // Write out cell
        write_cell_data(cell_idx, positions_start, velocities_start,positions_end, positions_end, all_status);
    } 

    //Delete data structures
    delete[] cell_particles;
    delete[] all_status;
    delete[] positions_start; 
    delete[] velocities_start;
    delete[] positions_end; 
    delete[] velocities_end;



}

void SimController :: setup_datawriter(){
    hid_t       file_id, dataset_id, dataspace_id, group;  /* identifiers */
    hsize_t  dims2D[4], dims3D[5];
    herr_t      status;
    char sidx[20], dset_name[30];

    dims2D[0] = sd->dim;
    dims2D[1] = sd->dim;
    dims2D[2] = sd->dim;
    dims2D[3] = N_particles;

    dims3D[0] = 3;
    dims3D[1] = sd->dim;
    dims3D[2] = sd->dim;
    dims3D[3] = sd->dim;
    dims3D[4] = N_particles;

    file_id = H5Fcreate(outname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dataspace_id = H5Screate_simple(4, dims2D, NULL);
    dataset_id = H5Dcreate2(file_id, "status", H5T_STD_I32BE,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    dataspace_id = H5Screate_simple(5, dims3D, NULL);
    dataset_id = H5Dcreate2(file_id, "position_start", H5T_NATIVE_FLOAT,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(5, dims3D, NULL);
    dataset_id = H5Dcreate2(file_id, "velocity_start", H5T_NATIVE_FLOAT,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(5, dims3D, NULL);
    dataset_id = H5Dcreate2(file_id, "position_end", H5T_NATIVE_FLOAT,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(5, dims3D, NULL);
    dataset_id = H5Dcreate2(file_id, "velocity_end", H5T_NATIVE_FLOAT,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    status = H5Fclose(file_id);
}

void SimController :: write_cell_data(int cell_idx, float * positions_start, float * velocities_start,float * positions_end, float * velocities_end, int * all_status){
    hid_t       file_id, dataset_id, dataspace_id, memspace_id;  /* identifiers */
    hsize_t  dims2D[4], dims3D[5], offset2D[4], offset3D[5];
    hsize_t stride[5] = {1,1,1,1,1}, block[5]={1,1,1,1,1};
    herr_t      status;
    int i, j, k;
    char sidx[20], dset_name[30];

    reverse_cellidx(cell_idx,sd->dim,i,j,k);
    

    dims2D[0] = 1;
    dims2D[1] = 1;
    dims2D[2] = 1;
    dims2D[3] = N_particles;

    dims3D[0] = 3;
    dims3D[1] = 1;
    dims3D[2] = 1;
    dims3D[3] = 1;
    dims3D[4] = N_particles;

    offset2D[0] = i;
    offset2D[1] = j;
    offset2D[2] = k;
    offset2D[3] = 0;

    offset3D[0] = 0;
    offset3D[1] = i;
    offset3D[2] = j;
    offset3D[3] = k;
    offset3D[4] = 0;



    // open existing file
    file_id = H5Fopen(outname, H5F_ACC_RDWR, H5P_DEFAULT);

    //create/close dataspace, dataset, write data, for status
    dataset_id = H5Dopen2(file_id, "status", H5P_DEFAULT);
    memspace_id = H5Screate_simple(4,dims2D, NULL); 
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D,
                                 stride, dims2D, block);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT,
                      all_status);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    //create/close dataspace, dataset, write data, for init location
    dataset_id = H5Dopen2(file_id, "position_start", H5P_DEFAULT);
    memspace_id = H5Screate_simple(5,dims3D, NULL); 
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D,
                                 stride, dims3D, block);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
                      positions_start);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    //create/close dataspace, dataset, write data, for init velocity
    dataset_id = H5Dopen2(file_id, "velocity_start", H5P_DEFAULT);
    memspace_id = H5Screate_simple(5,dims3D, NULL); 
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D,
                                 stride, dims3D, block);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
                      velocities_start);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    //create/close dataspace, dataset, write data, for init location
    dataset_id = H5Dopen2(file_id, "position_end", H5P_DEFAULT);
    memspace_id = H5Screate_simple(5,dims3D, NULL); 
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D,
                                 stride, dims3D, block);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
                      positions_end);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);


    //create/close dataspace, dataset, write data, for init velocity
    dataset_id = H5Dopen2(file_id, "velocity_end", H5P_DEFAULT);
    memspace_id = H5Screate_simple(5,dims3D, NULL); 
    dataspace_id = H5Dget_space (dataset_id);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D,
                                 stride, dims3D, block);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT,
                      velocities_end);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    //close group and file
    status = H5Fclose(file_id);
}

bool SimController :: eval_cell(int cell_idx){
    bool evaluate = false;
    float x,y,z,r;
    int i,j,k;

    reverse_cellidx(cell_idx, sd->dim, i,j,k);

    x = sd->x[i];
    y = sd->y[j];
    z = sd->z[k];
    r = pow(x*x+y*y+z*z,0.5);

    //want < 2.5 RM, > 1 RM
    if ((r > 3390+300-100) && (r < 3390+800+100)){evaluate=true;}
    //if ((r > 3390+800-100) && (r < 3390+800+100)){evaluate=true;}
    
    //if ((i == 59)&&(j==43)&&(k==59)){evaluate=true;}
    //else{evaluate=false;}
    
    /*
    evaluate=false;
    if (i == 59){
        if (j == 59){
            if ((k>=75) && (k<=82)){ evaluate=true;}
            if ((k<=43) && (k>=37)){ evaluate=true;}
        }
        if (k == 59){
            if ((j>=75) && (j<=82)){ evaluate=true;}
            if ((j<=43) && (j>=37)){ evaluate=true;}
        }
    }
    if ((j == 59) && (k == 59)){
        if ((i>=75)&&(i<=82)){evaluate=true;}
    }
    */
    return evaluate;
}

void SimController :: set_particle_pop(float mass, float charge, float temperature){
    pop_m = mass;
    pop_q = charge;
    pop_T = temperature;

    pop_set=true;
}
