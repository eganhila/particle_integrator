#define const_Bz
#include <iostream>
#include "H5Cpp.h"
using namespace H5;

int read_simulation_data(){
	const H5std_string FILE_NAME("/Users/hilaryegan/Data/MagneticField/PrelimAllEnd/B_0nT.h5");
	H5std_string;// DATASET_NAME()
	const int N=120;
	float magnetic_field[3,N,N,N];
	float coordinates[3,N];
	hsize_t  offset_out[3];   // hyperslab offset in memory
    hsize_t  count_out[3]; 
	
	hsize_t     dimsm[3] = {N,N,N};              /* memory space dimensions */
    DataSpace memspace( RANK_OUT, dimsm );


	H5File file( FILE_NAME, H5F_ACC_RDONLY );


	DataSet dataset = file.openDataSet("magnetic_field_x");
	DataSpace dataspace = dataset.getSpace();
	memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out )
	dataset.read(magnetic_field[0], PredType::NATIVE_FLOAT, memspace, dataspace );

	DataSet dataset = file.openDataSet("magnetic_field_y");
	DataSpace dataspace = dataset.getSpace();
	memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out )
	dataset.read(magnetic_field[1], PredType::NATIVE_FLOAT, memspace, dataspace );


	DataSet dataset = file.openDataSet("magnetic_field_z");
	DataSpace dataspace = dataset.getSpace();
	memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out )
	dataset.read(magnetic_field[2], PredType::NATIVE_FLOAT, memspace, dataspace );
}


struct Particle{
	float state[6];
	float mass;     // velocity
};

struct Derivative {
	float d[6]={0};
};


float acceleration( const Particle & particle, double t, int axis){

    #ifdef spring
        const float k = 15.0f;
    	const float b = 0.1f;

     	return -k * particle.state[axis] - b * particle.state[3+axis];
    #endif

     #ifdef const_Bz
     	const float B[3] = {0,0,5};
     	return -1*particle.state[3+(axis+1)%3]*B[(axis+2)%3]+particle.state[3+(axis+2)%3]*B[(axis+1)%3];
     #endif

}

Derivative evaluate( const Particle & particle_init, 
                     double t, 
                     float dt, 
                     const Derivative & d ){
    Particle temp;

    for (int i =0; i<6; i++){
    	temp.state[i] = particle_init.state[i] + d.d[i] * dt;
	}

    Derivative output;
    for (int i =0; i<3; i++){
    	output.d[i] = temp.state[3+i];
    	output.d[3+i] = acceleration( temp, t+dt, i);
	}
    return output;
}


void integrate( Particle & particle, 
                double t, 
                float dt ){

    Derivative a,b,c,d;//,k;
    //for (int i =0; i<3; i++){
    //	k.d[i] = particle.state[3+i];
    //	k.d[3+i] = acceleration(particle, t, i);
	//}

    a = evaluate( particle, t, 0.0f, Derivative());
    b = evaluate( particle, t, dt*0.5f, a );
    c = evaluate( particle, t, dt*0.5f, b );
    d = evaluate( particle, t, dt, c );

    float ddt[6];

    for (int i =0; i<6; i++){
    	ddt[i] = 1.0f/6.0f* (a.d[i]+2.0f* (b.d[i]+c.d[i])+d.d[i]);
    	particle.state[i] = particle.state[i] + ddt[i] * dt;
	}
}


int main(){
	Particle particle;
	float t=0, dt=0.001, T=2;

	particle.state[0] = 0;
	particle.state[1] = 0;
	particle.state[2] = 0;
	particle.state[3] = 1;
	particle.state[4] = 1;
	particle.state[5] = 1;

	read_simulation_data()

	while (t < T){
		integrate(particle, t, dt);
		std::cout << particle.state[0]<<" " << particle.state[1]<<"\n";
		t += dt;
	}
	


	return 0;
}