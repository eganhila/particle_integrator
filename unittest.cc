#ifdef GTEST
#include "gtest/gtest.h"
#include "sim_dat.h"
#include "interpolate.h"

TEST(SimDatTest, Constructor) {
    SimDat sd(10);
    
    EXPECT_EQ(10, sd.dim);
    EXPECT_EQ(10, sizeof(sd.x)/sizeof(float));
    EXPECT_EQ(10*10*10, sizeof(sd.Bx)/sizeof(float));
}

GTEST_API_ int main(int argc, char **argv) {
  printf("Running main() from %s\n", __FILE__);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#endif
