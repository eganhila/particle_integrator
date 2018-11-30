#ifdef GTEST
#include "gtest/gtest.h"
#include "sim_dat.h"
#include "interpolate.h"
#include "integrator.h"

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
//------------------------ INTEGRATOR TESTS -----------------------------///

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
