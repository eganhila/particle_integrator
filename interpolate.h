#include "sim_dat.h"

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

void LinearInterpolate(const float * c0, const float * c1, float xd, float * c);
bool TrilinearInterpolate(const float * pos, const SimDat & sd, float * pss);
bool getCellIdx(const SimDat & sd, const float * pos, int *idx);

#endif /* INTERPOLATE_H */
