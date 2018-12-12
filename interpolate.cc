#include "interpolate.h"
#include "sim_dat.h"
#include <iostream>
void LinearInterpolate(const float * c0, const float * c1, float xd, float * c){
    for (int i=0; i<6; i++){
        c[i] = c0[i]*(1-xd)+c1[i]*xd;
    }
}

bool getCellIdx(const SimDat & sd, const float * pos, int *idx){
    // Find cell we are in
    int i0=0, j0=0, k0=0; //leftmost corner

    while (i0<sd.dim){ 
        if (pos[0]< sd.x[i0]){break;}
        i0+=1;}
    while (j0<sd.dim){ 
        if (pos[1]< sd.y[j0]){break;}
        j0+=1;}
    while (k0<sd.dim){ 
        if (pos[2]< sd.z[k0]){break;}
        k0+=1;}

    idx[0] = i0-1;
    idx[1] = j0-1;
    idx[2] = k0-1;

    if ((i0>=sd.dim) || (j0>=sd.dim) || k0>=sd.dim){return false;}

    return true;
}

bool TrilinearInterpolate(const float * pos, const SimDat & sd, float * pss){

    int idx[3];

    if (!getCellIdx(sd, pos, idx)){return false;} 


    //Find distance to cell
    float xd,yd,zd;
    xd = (pos[0]-sd.x[idx[0]])/(sd.x[idx[0]+1]-sd.x[idx[0]]);
    yd = (pos[1]-sd.y[idx[1]])/(sd.y[idx[1]+1]-sd.y[idx[1]]);
    zd = (pos[2]-sd.z[idx[2]])/(sd.z[idx[2]+1]-sd.z[idx[2]]);

    int i0=idx[0], j0=idx[1], k0=idx[2];



    // Find Data of neighbor cells
    float c000[6], c001[6], c010[6], c011[6], c100[6], c101[6], c110[6], c111[6];
    
    sd.GetSimState(i0  , j0  , k0  , c000);
    sd.GetSimState(i0  , j0  , k0+1, c001);
    sd.GetSimState(i0  , j0+1, k0  , c010);
    sd.GetSimState(i0  , j0+1, k0+1, c011);
    sd.GetSimState(i0+1, j0  , k0  , c100);
    sd.GetSimState(i0+1, j0  , k0+1, c101);
    sd.GetSimState(i0+1, j0+1, k0  , c110);
    sd.GetSimState(i0+1, j0+1, k0+1, c111);
    

    // Interpolate data
    float c00[6], c01[6], c10[6], c11[6], c0[6], c1[6];

    LinearInterpolate(c000, c100, xd, c00);
    LinearInterpolate(c001, c101, xd, c01);
    LinearInterpolate(c010, c110, xd, c10);
    LinearInterpolate(c011, c111, xd, c11);

    LinearInterpolate(c00, c10, yd, c0);
    LinearInterpolate(c01, c11, yd, c1);

    LinearInterpolate(c0, c1, zd, pss);

    return true;
}

