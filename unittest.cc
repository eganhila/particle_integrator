#ifdef GTEST
#include "gtest/gtest.h"
#include "sim_dat.h"
#include "interpolate.h"
#include "integrator.h"
#include "acceleration.h"

void setupFakeSD(SimDat & sd){
    int i,j,k,d=sd.dim, idx;

    for (i=0;i<d;i++){
        for (j=0;j<d;j++){
            for (k=0;k<d;k++){
                idx = i*d*d+j*d+k;
                sd.Bx[idx] = i;
                sd.By[idx] = j;
                sd.Bz[idx] = k;
                sd.Ex[idx] = 2;
                sd.Ey[idx] = 1;
                sd.Ez[idx] = 0;
            }}
        sd.x[i] = i;
        sd.y[i] = i/2.0;
        sd.z[i] = i+1;
    }
}
void setupZeroEFakeSD(SimDat & sd){
    int i,j,k,d=sd.dim, idx;

    for (i=0;i<d;i++){
        for (j=0;j<d;j++){
            for (k=0;k<d;k++){
                idx = i*d*d+j*d+k;
                sd.Bx[idx] = i;
                sd.By[idx] = j;
                sd.Bz[idx] = k;
                sd.Ex[idx] = 0;
                sd.Ey[idx] = 0;
                sd.Ez[idx] = 0;
            }}
        sd.x[i] = i;
        sd.y[i] = i/2.0;
        sd.z[i] = i+1;
    }
}

void setupConstSD(SimDat & sd){
    int d = sd.dim;
    int d3 = d*d*d;
        for (int i=0; i<(d3); i++){
            sd.Bx[i] = 0;
            sd.By[i] = 0;
            sd.Bz[i] = 2;
            sd.Ex[i] = 0;
            sd.Ey[i] = 0;
            sd.Ez[i] = 0;
        }
        for (int i=0; i<(d); i++){
            sd.x[i] = -4+8.0*i/d;
            sd.y[i] = -4+8.0*i/d;
            sd.x[i] = -4+8.0*i/d;
        }
        sd.set_bounds();
}
//------------------------ ACCELERATION TESTS-----------------------------///

TEST(AccelerationTests, ConstBz_Zero){
    Particle p;
    float acc[3];
    SimDat sd(3);

    for (int i=0; i<6; i++){p.state[i]= 0;}
    
    constBz_accel(p, sd, acc);

    EXPECT_EQ(acc[0], 0.0);
    EXPECT_EQ(acc[1], 0.0);
    EXPECT_EQ(acc[2], 0.0);
}

TEST(AccelerationTests, ConstBz_NonZero){
    Particle p;
    float acc[3];
    SimDat sd(3);

    for (int i=0; i<6; i++){p.state[i]= 0;}
    p.state[3] = 1;
    p.state[5] = 3;
    
    constBz_accel(p, sd, acc);

    EXPECT_EQ(acc[0], 0.0);
    EXPECT_EQ(acc[1], 2);
    EXPECT_EQ(acc[2], 0.0);
}

TEST(AccelerationTests, SimDat_Zero){
    Particle p;
    float acc[3];
    SimDat sd(3);
    setupZeroEFakeSD(sd);

    for (int i=0; i<3; i++){
        p.state[i]= 1;
        p.state[i+3] = 0;}

    simDat_accel(p, sd, acc);

    EXPECT_EQ(acc[0], 0.0);
    EXPECT_EQ(acc[1], 0.0);
    EXPECT_EQ(acc[2], 0.0);


}

TEST(AccelerationTests, SimDat_ConstBz){
    Particle p;
    float acc[3];
    SimDat sd(3);
    int d = 3;

    for (int i=0; i<3; i++){p.state[i]= 1;}
    p.state[3] = 1;
    p.state[4] = 0;
    p.state[5] = 3;

    setupFakeSD(sd);
    for (int i=0; i<27; i++){
        sd.Bx[i] = 0;
        sd.By[i] = 0;
        sd.Bz[i] = 2;
        sd.Ex[i] = 0;
        sd.Ey[i] = 0;
        sd.Ez[i] = 0;
    }

    simDat_accel(p, sd, acc);

    EXPECT_EQ(acc[0], 0.0);
    EXPECT_EQ(acc[1], 2);
    EXPECT_EQ(acc[2], 0.0);

}



//------------------------ INTEGRATOR TESTS -----------------------------///

TEST(IntegratorTest, IntegrateStationary){
    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = 0;}
    Integrator intg(p, 1,5, 1);
    intg.set_accel(constBz_accel);

    intg.integrate();

    EXPECT_FLOAT_EQ(intg.t_final, 5);
    EXPECT_FLOAT_EQ(p.state[0], 0);
    EXPECT_FLOAT_EQ(p.state[1], 0);
    EXPECT_FLOAT_EQ(p.state[2], 0);
    EXPECT_FLOAT_EQ(p.state[3], 0);
    EXPECT_FLOAT_EQ(p.state[4], 0);
    EXPECT_FLOAT_EQ(p.state[5], 0);

}

TEST(IntegratorTest, IntegrateLinear){
    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = 1;}
    Integrator intg(p, 1,5, 2);
    intg.set_accel(zero_accel);

    intg.integrate();

    EXPECT_FLOAT_EQ(intg.t_final, 5);
    EXPECT_FLOAT_EQ(p.state[0], 5);
    EXPECT_FLOAT_EQ(p.state[1], 5);
    EXPECT_FLOAT_EQ(p.state[2], 5);
    EXPECT_FLOAT_EQ(p.state[3], 1);
    EXPECT_FLOAT_EQ(p.state[4], 1);
    EXPECT_FLOAT_EQ(p.state[5], 1);
}

TEST(IntegratorTest, IntegrateConstBzAcc){
    Particle p;
    p.state[0] = 0;
    p.state[1] = 0;
    p.state[2] = 0;
    p.state[3] = 0;
    p.state[4] = 1;
    p.state[5] = 0;
    Integrator intg(p, 1,5, 0.01);
    intg.set_accel(constBz_accel);

    intg.integrate();

    EXPECT_FLOAT_EQ(p.state[0], -0.57274985);
    EXPECT_FLOAT_EQ(p.state[1], 0.49467942);
    EXPECT_FLOAT_EQ(p.state[2], 0);
    EXPECT_FLOAT_EQ(p.state[3], -0.98935884);
    EXPECT_FLOAT_EQ(p.state[4], -0.14550011);
    EXPECT_FLOAT_EQ(p.state[5], 0);
}

TEST(IntegratorTest, IntegrateStationarySimdat){
    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = 0;}
    Integrator intg(p, 1,5, 1);
    intg.set_accel(simDat_accel);
    SimDat sd(3);
    setupZeroEFakeSD(sd);
    intg.set_sd(sd);

    intg.integrate();

    //having trouble comparing floats very near zero, so added 1
    EXPECT_FLOAT_EQ(intg.t_final, 5);
    EXPECT_FLOAT_EQ(p.state[0]+1,1);
    EXPECT_FLOAT_EQ(p.state[1]+1,1);
    EXPECT_FLOAT_EQ(p.state[2]+1,1);
    EXPECT_FLOAT_EQ(p.state[3]+1,1);
    EXPECT_FLOAT_EQ(p.state[4]+1,1);
    EXPECT_FLOAT_EQ(p.state[5]+1,1);

}

TEST(IntegratorTest, IntegrateConstBzAccSimDat){
    Particle p;
    p.state[0] = 0;
    p.state[1] = 0;
    p.state[2] = 0;
    p.state[3] = 0;
    p.state[4] = 1;
    p.state[5] = 0;
    Integrator intg(p, 1,5, 0.01);
    SimDat sd(20);
    setupConstSD(sd);
    intg.set_sd(sd);

    intg.set_accel(simDat_accel);
    intg.integrate();

    EXPECT_FLOAT_EQ(p.state[0], -0.57274985);
    EXPECT_FLOAT_EQ(p.state[1], 0.49467942);
    EXPECT_FLOAT_EQ(p.state[2], 0);
    EXPECT_FLOAT_EQ(p.state[3], -0.98935884);
    EXPECT_FLOAT_EQ(p.state[4], -0.14550011);
    EXPECT_FLOAT_EQ(p.state[5], 0);
}


TEST(IntegratorTest, EvaluateDerivZero){
    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = 0;}
    Integrator intg(p, 1,5, 2);
    intg.set_accel(constBz_accel);

    float d_in[6]={0}, d_out[6];

    intg.evaluate_derivative(2, 1, d_in, d_out);

    EXPECT_EQ(d_out[0], 0);
    EXPECT_EQ(d_out[1], 0);
    EXPECT_EQ(d_out[2], 0);
    EXPECT_EQ(d_out[3], 0);
    EXPECT_EQ(d_out[4], 0);
    EXPECT_EQ(d_out[5], 0);

}
TEST(IntegratorTest, EvaluateDerivBase){
    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = 0;}
    p.state[3] = 1;
    p.state[5] = 3;
    Integrator intg(p, 1,6,  1);
    intg.set_accel(constBz_accel);

    float d_in[6]={0}, d_out[6];

    intg.evaluate_derivative(2, 1, d_in, d_out);

    EXPECT_EQ(d_out[0], 1);
    EXPECT_EQ(d_out[1], 0);
    EXPECT_EQ(d_out[2], 3);
    EXPECT_EQ(d_out[3], 0);
    EXPECT_EQ(d_out[4], 2);
    EXPECT_EQ(d_out[5], 0);

}
TEST(IntegratorTest, EvaluateDerivSmallDt){
    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = 0;}
    p.state[3] = 1;
    p.state[5] = 3;
    Integrator intg(p, 1,5, 0.5);
    intg.set_accel(constBz_accel);

    float d_in[6]={0}, d_out[6];

    intg.evaluate_derivative(2, 0.5, d_in, d_out);

    EXPECT_EQ(d_out[0], 0.5);
    EXPECT_EQ(d_out[1], 0);
    EXPECT_EQ(d_out[2], 1.5);
    EXPECT_EQ(d_out[3], 0);
    EXPECT_EQ(d_out[4], 1);
    EXPECT_EQ(d_out[5], 0);

}


TEST(IntegratorTest, EvaluateDerivDin){

    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = 0;}
    p.state[3] = 1;
    p.state[5] = 3;
    Integrator intg(p, 1,5, 0.5);
    intg.set_accel(constBz_accel);

    float d_in[6]={1,0,0,1,0,0}, d_out[6];

    intg.evaluate_derivative(2, 1, d_in, d_out);

    EXPECT_EQ(d_out[0], 2);
    EXPECT_EQ(d_out[1], 0);
    EXPECT_EQ(d_out[2], 3);
    EXPECT_EQ(d_out[3], 0);
    EXPECT_EQ(d_out[4], 4);
    EXPECT_EQ(d_out[5], 0);
}

TEST(IntegratorTest, RK4MethodLinear){

    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = 0;}
    p.state[3] = 1;
    p.state[5] = 3;
    Integrator intg(p, 0,5, 1);
    intg.set_accel(constBz_accel);

    intg.integrate_step();


    EXPECT_FLOAT_EQ(p.state[0], 1.0f/3.0f);
    EXPECT_FLOAT_EQ(p.state[1], 2.0f/3.0f);
    EXPECT_FLOAT_EQ(p.state[2], 3);
    EXPECT_FLOAT_EQ(p.state[3], -1/3.0f);
    EXPECT_FLOAT_EQ(p.state[4], 2.0f/3.0f);
    EXPECT_FLOAT_EQ(p.state[5], 3.0);
    EXPECT_FLOAT_EQ(intg.t, 1);

}
TEST(IntegratorTest, RK4MethodLinearSmallDt){

    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = 0;}
    p.state[3] = 1;
    p.state[5] = 3;
    Integrator intg(p, 0,4,  0.5);
    intg.set_accel(constBz_accel);

    intg.integrate_step();


    EXPECT_FLOAT_EQ(p.state[0], 5.0f/12.0f);
    EXPECT_FLOAT_EQ(p.state[1], 11.0f/48.0f);
    EXPECT_FLOAT_EQ(p.state[2], 1.5);
    EXPECT_FLOAT_EQ(p.state[3], 13.0f/24.0);
    EXPECT_FLOAT_EQ(p.state[4], 5.0f/6.0f);
    EXPECT_FLOAT_EQ(p.state[5], 3.0);

    EXPECT_EQ(intg.t, 0.5);
}


TEST(IntegratorTest, ConstructorIntegers){
    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = i;}
    Integrator intg(p, 1,4,  2);

    EXPECT_EQ(intg.t0, 1);
    EXPECT_EQ(intg.dt, 2);
    EXPECT_EQ(intg.t, 1);
    EXPECT_EQ(intg.t_final, 4);
    EXPECT_EQ(intg.particle->state[0], 0);
    EXPECT_EQ(intg.particle->state[4], 4);
    EXPECT_EQ(intg.has_sd, false);

}

TEST(IntegratorTest, ConstructorNonIntegers){
    Particle p;
    for (int i=0; i<6; i++){ p.state[i] = float(i)/2.0;}
    Integrator intg(p, 0.5,50.6,  0.01);

    EXPECT_EQ(intg.t0, 0.5);
    EXPECT_EQ(intg.dt, float(0.01));
    EXPECT_EQ(intg.t, 0.5);
    EXPECT_FLOAT_EQ(intg.t_final, 50.6);
    EXPECT_EQ(intg.particle->state[0], 0);
    EXPECT_EQ(intg.particle->state[3], 1.5);
}

TEST(IntegratorTest, SetSD){
    Particle p;
    Integrator intg(p, 0.5, 1,  0.01);
    SimDat sd(3);

    intg.set_sd(sd);

    EXPECT_EQ(intg.sd, &sd);
    EXPECT_EQ(intg.has_sd, true);
}

TEST(IntegratorTest, SetAccel){
    Particle p;
    Integrator intg(p, 0.5,5, 0.01);

    intg.set_accel(simDat_accel);

    EXPECT_EQ(intg.acc_func, &simDat_accel);

}




///----------------------- INTERPOLATE TESTS ----------------------------///

TEST(TrilinearInterpolateTest, Normal){
    SimDat sd(3);
    setupFakeSD(sd);

    float pos[3] = {0.5,0.75,2.5};
    float pss_out[6];

    TrilinearInterpolate(pos, sd, pss_out);
    EXPECT_EQ(pss_out[0], 0.5);
    EXPECT_EQ(pss_out[1], 1.5);
    EXPECT_EQ(pss_out[2], 1.5);
    EXPECT_EQ(pss_out[3], 2);
    EXPECT_EQ(pss_out[4], 1);
    EXPECT_EQ(pss_out[5], 0);

}

TEST(GetCellIdxTest, Normal){
    SimDat sd(3);
    setupFakeSD(sd);

    float pos[3] = {0.5,0.75,3.5};
    int idx[3];
    getCellIdx(sd, pos, idx);

    EXPECT_EQ(idx[0], 0);
    EXPECT_EQ(idx[1], 1);
    EXPECT_EQ(idx[2], 2);
}

TEST(GetCellIdxTest, OnBound){
    SimDat sd(3);
    setupFakeSD(sd);

    float pos[3] = {1,0,3};
    int idx[3];
    getCellIdx(sd, pos, idx);

    EXPECT_EQ(idx[0], 1);
    EXPECT_EQ(idx[1], 0);
    EXPECT_EQ(idx[2], 2);
}

// Stubs b/c hopefully we catch this outside
TEST(GetCellIdxTest, OutOfBoundsLT){
    
}

// Stubs b/c hopefully we catch this outside
TEST(GetCellIdxTest, OutOfBoundsGT){
    
}

TEST(LinearInterpolateTest, Normal){
    float c0[6], c1[6], c[6];
    for (int i=0; i<6; i++){
        c0[i] = 0;
        c1[i] = i;}

    LinearInterpolate(c0,c1,0.5,c);
    EXPECT_EQ(c[1], 0.5);
    EXPECT_EQ(c[2], 1);

    LinearInterpolate(c0,c1,0.75,c);
    EXPECT_EQ(c[3], 2.25);
    EXPECT_EQ(c[4], 3);

    LinearInterpolate(c0,c1,1,c);
    EXPECT_EQ(c[0], 0);
    EXPECT_EQ(c[5], 5);
}






////-------------- SIM DAT TESTS ----------------------------------------///
TEST(SimDatTest, Constructor) {
    SimDat sd(10);
    EXPECT_EQ(10, sd.dim);
    //Can't test sizeof new arrays, not stored b/c dynamic
}

TEST(SimDatTest, GetSimState) {
    SimDat sd(3);
    setupFakeSD(sd);
    
    float out[6];

    sd.GetSimState(0,0,0,out);
    EXPECT_EQ(out[0],0);
    EXPECT_EQ(out[1],0);
    EXPECT_EQ(out[2],0);
    EXPECT_EQ(out[3],2);
    
    sd.GetSimState(2,2,2,out);
    EXPECT_EQ(out[0],2);
    EXPECT_EQ(out[1],2);
    EXPECT_EQ(out[2],2);
    EXPECT_EQ(out[5],0);

    sd.GetSimState(0,1,2,out);
    EXPECT_EQ(out[0],0);
    EXPECT_EQ(out[1],1);
    EXPECT_EQ(out[2],2);
    EXPECT_EQ(out[4],1);
}

TEST(ReadSDTest, coordinates){
    SimDat sd(3);
    read_simulation_data(sd, "test/read_test.h5");

    EXPECT_EQ(sd.x[0], 0);
    EXPECT_EQ(sd.x[1], 1);
    EXPECT_EQ(sd.x[2], 2);

    EXPECT_EQ(sd.y[0], 0);
    EXPECT_EQ(sd.y[1], 0.5);
    EXPECT_EQ(sd.y[2], 1);
    
    EXPECT_EQ(sd.z[0], 1);
    EXPECT_EQ(sd.z[1], 2);
    EXPECT_EQ(sd.z[2], 3);
}

TEST(ReadSDTest, data){
    SimDat sd(3);
    read_simulation_data(sd, "test/read_test.h5");
    float out[6];

    sd.GetSimState(0,0,0,out);
    EXPECT_EQ(out[0],0);
    EXPECT_EQ(out[1],0);
    EXPECT_EQ(out[2],0);
    EXPECT_EQ(out[3],2);
    
    sd.GetSimState(2,2,2,out);
    EXPECT_EQ(out[0],2);
    EXPECT_EQ(out[1],2);
    EXPECT_EQ(out[2],2);
    EXPECT_EQ(out[5],0);

    sd.GetSimState(0,1,2,out);
    EXPECT_EQ(out[0],0);
    EXPECT_EQ(out[1],1);
    EXPECT_EQ(out[2],2);
    EXPECT_EQ(out[4],1);
    
}

GTEST_API_ int main(int argc, char **argv) {
  printf("Running main() from %s\n", __FILE__);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#endif
