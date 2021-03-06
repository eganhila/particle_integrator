#include "sim_dat.h"
#include <string>
#include "hdf5.h"
#define FILE            "/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_0nT.h5"

int read_simulation_data(SimDat& sd, const char* file_name){
    hid_t       file, space, dset, dsetx, dsety, dsetz;          /* Handles */
    herr_t      status;
    int            i, j, k;

    file = H5Fopen (file_name, H5F_ACC_RDONLY, H5P_DEFAULT);

    //x,y,z
    dsetx = H5Dopen (file, "x", H5P_DEFAULT);
    dsety = H5Dopen (file, "y", H5P_DEFAULT);
    dsetz = H5Dopen (file, "z", H5P_DEFAULT);

    status = H5Dread(dsetx, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.x);
    status = H5Dread(dsety, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.y);
    status = H5Dread(dsetz, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.z);

    status = H5Dclose (dsetx);
    status = H5Dclose (dsety);
    status = H5Dclose (dsetz);
    
    //mag field
    dsetx = H5Dopen (file, "magnetic_field_x", H5P_DEFAULT);
    dsety = H5Dopen (file, "magnetic_field_y", H5P_DEFAULT);
    dsetz = H5Dopen (file, "magnetic_field_z", H5P_DEFAULT);

    status = H5Dread(dsetx, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Bx);
    status = H5Dread(dsety, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.By);
    status = H5Dread(dsetz, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Bz);

    status = H5Dclose (dsetx);
    status = H5Dclose (dsety);
    status = H5Dclose (dsetz);
    
    dsetx = H5Dopen (file, "motional_electric_field_x", H5P_DEFAULT);
    dsety = H5Dopen (file, "motional_electric_field_y", H5P_DEFAULT);
    dsetz = H5Dopen (file, "motional_electric_field_z", H5P_DEFAULT);

    status = H5Dread(dsetx, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Ex);
    status = H5Dread(dsety, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Ey);
    status = H5Dread(dsetz, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sd.Ez);

    status = H5Dclose (dsetx);
    status = H5Dclose (dsety);
    status = H5Dclose (dsetz);

    status = H5Fclose (file);

    sd.set_bounds();

    return 0;
}

void SimDat::set_bounds(){
    bbox[0] = x[0];
    bbox[1] = y[0];
    bbox[2] = z[0];
    bbox[3] = x[dim-1];
    bbox[4] = y[dim-1];
    bbox[5] = z[dim-1];

    dx=(x[dim-1]-x[0])/dim;
}
