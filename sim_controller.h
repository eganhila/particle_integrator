#ifndef SIM_C
#define SIM_C

#include "sim_dat.h"
#include "integrator.h"

class SimController {
    public:
        int N_particles;
        SimDat * sd = NULL;
        float t_final, dt;
        bool created_particles=false;
        Particle * init_pop;

    void setup();
    void draw_init_pop();
    void run();
    void write_out();


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
