#include "sim_dat.h"
#include "integrator.h"
#include "sim_controller.h"

void SimController :: setup(){
    
    init_pop = new Particle[N_particles];
    draw_init_pop();
}

void SimController :: draw_init_pop(){
}

void SimController :: run(){
    write_out();
}

void SimController :: write_out(){
}
