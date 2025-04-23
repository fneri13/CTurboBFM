#include <iostream>
#include "../../../include/CFluid.hpp"
#include "gtest/gtest.h"


TEST(CFluidTest, TestStaticEnergy_p_rho) {
    FloatType gamma {2.0}, R {1.0};
    CFluid fluid(gamma, R);

    std::vector<FloatType> testDensity {1.0};
    std::vector<FloatType> testPressure {1.0};
    std::vector<FloatType> expectedStaticEnergy {1.0};

    for (int i=0; i<testDensity.size(); i++){
        ASSERT_DOUBLE_EQ(fluid.computeStaticEnergy_p_rho(testPressure[i], testDensity[i]), expectedStaticEnergy[i]);
    }
    
}

TEST(CFluidTest, TestSoundSpeed_p_rho) {
    FloatType gamma {2.0}, R {1.0};
    CFluid fluid(gamma, R);
    
    FloatType rho {1.0};
    FloatType p {1.0};

    FloatType expectedValue {std::sqrt(2.0)};

    std::vector<FloatType> computedStaticEnergy;
    ASSERT_DOUBLE_EQ(fluid.computeSoundSpeed_p_rho(p, rho), expectedValue);
    
}



int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}