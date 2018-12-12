#include "sim_dat.h"
#include "integrator.h"
#include "sim_controller.h"
#include "acceleration.h"
#include <random>



Particle SimController :: draw_particle(){
    //Chapman initial pop

    /*float x, y, z, H, a0, a1, a, R0;
    float noon_f=1.0, night_f=0.1;

    //setup randoming
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> uniform_dist(0, 1);

    // 2.213 km/K assuming all mars things
    H = 2.213*pop_t/pop_m;


    x = uniform_dist(e2);
    y = uniform_dist(e2);
    z = uniform_dist(e2);
    */

    Particle p;
    return p;
}

void SimController :: run(){
    int cell_idx;
    int p_idx;
    Particle * cell_particles;
    int * all_status;

    // Setup data structures that we will
    // keep re-using
    cell_particles = new Particle[N_particles]; 
    all_status = new int[N_particles];

    // Iterate over all cells
    for (cell_idx=0; cell_idx<sd->dim3; cell_idx++){

        //Evaluate if we want to check this cell
        if (!eval_cell(cell_idx)){continue;}
        

        // Each cell needs N_particles, this is what we're going
        // to parallelize over
        for (p_idx=0; p_idx< N_particles; p_idx++){

            //draw particle from distribution
            cell_particles[p_idx] = draw_particle();

            //Setup integrator
            Integrator intg(cell_particles[p_idx], 1,100, 0.1);
            intg.set_sd(*sd);
            intg.set_accel(simDat_accel);
            
            //Particle integrated (work done here)
            all_status[p_idx] = intg.integrate();
            
        }

        // Write out cell
        write_cell_data(cell_particles, all_status);
    } 

    //Delete data structures
    delete[] cell_particles;
    delete[] all_status;



}

void SimController :: write_cell_data(Particle * cell_particles, int * all_status){
}

bool SimController :: eval_cell(int cell_idx){

    //For now we're just going to evaluate a single cell
    if (cell_idx == 0){ return true;}
    
    return false;
}

void SimController :: set_particle_pop(float mass, float charge, float temperature){
    pop_m = mass;
    pop_q = charge;
    pop_T = temperature;

    pop_set=true;
}
