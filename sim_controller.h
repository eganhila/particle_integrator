#ifndef SIM_C
#define SIM_C

#include "sim_dat.h"
#include "integrator.h"
#define MAX_STEPS 100000

class SimController {
    private:
        float pop_m, pop_q, pop_T;
        SimDat * sd = NULL;
        float t_final, dt;
        bool created_particles=false, pop_set=false;
    public:
        int N_particles;
        Particle * init_pop;
        const char * outname;
        bool uniform_E = false;

    Particle draw_particle(int cell_idx);
    void run_sim();
    void run_cell(int i, int j, int k);
    void write_cell_data(int cell_idx, float * positions, float * velocities, int * all_status);
    void set_particle_pop(float mass, float charge, float temperature);
    bool eval_cell(int cell_idx);
    void setup_datawriter();
    void setup_particlewriter();
    void write_particle(int p_idx, int N, float  pdata[6][MAX_STEPS]);


    SimController(int new_N_particles, 
                  float new_t_final,
                  float new_dt,
                  SimDat & new_sd,
                  const char * new_outname) :
        N_particles(new_N_particles), t_final(new_t_final),
        dt(new_dt), sd(&new_sd), outname(new_outname){}

    ~SimController(){
        if (created_particles){
            delete[] init_pop;
        }
    }
};

#endif
