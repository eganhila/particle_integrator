
struct Particle{
        float state[6];
        float mass;     
};

struct Derivative {
        float d[6]={0};
};

float acceleration(const Particle & particle, double t, int axis);
Derivative evaluate(const Particle &particle_init, double t, float dt, const Derivative & d);
void integrate(Particle & particle, double t, float dt);
