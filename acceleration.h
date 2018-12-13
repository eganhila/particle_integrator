#ifndef ACCEL_H
#define ACCEL_H

#include <iostream>
#include "integrator.h"
#include "sim_dat.h"
#include "interpolate.h"



bool zero_accel(const Particle & particle, const SimDat &sd, float * acc);
bool constBz_accel(const Particle & particle, const SimDat &sd, float * acc);
bool simDat_accel(const Particle & particle, const  SimDat &sd, float * acc);

#endif
