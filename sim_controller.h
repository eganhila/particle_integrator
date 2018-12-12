#ifndef SIM_C
#define SIM_C

#include "sim_dat.h"
#include "integrator.h"

class SimController {
    private:
        float pop_m, pop_q, pop_T;
        SimDat * sd = NULL;
        float t_final, dt;
        bool created_particles=false, pop_set=false;
    public:
        int N_particles;
        Particle * init_pop;

    Particle draw_particle();
    void run();
    void write_cell_data(int cell_idx, float * positions, float * velocities, int * all_status);
    void set_particle_pop(float mass, float charge, float temperature);
    bool eval_cell(int cell_idx);


    SimController(int new_N_particles, 
                  float new_t_final,
                  float new_dt,
                  SimDat & new_sd) :
        N_particles(new_N_particles), t_final(new_t_final),
        dt(new_dt), sd(&new_sd){}

    ~SimController(){
        if (created_particles){
            delete[] init_pop;
        }
    }
};

#endif
