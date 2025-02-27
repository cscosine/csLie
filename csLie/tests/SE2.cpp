#include <catch2/catch_test_macros.hpp>

#include <csLie/SE2.h>
#include <csLie/SE2Alg.h>
#include <csLie/SE2Exp.h>
#include <csLie/SE2Params.h>

static const int N_RANDOM = 1000;

TEST_CASE( "SE2 params2Pose2params", "[SE2]" ) {
    
  for (int i = 0; i < N_RANDOM; i++) {
    csLie::SE2Paramsd par;
    par().setRandom();
    auto pose = csLie::params2Pose(par);
    auto par2 = csLie::pose2Params(pose);
    REQUIRE((par2() - par()).squaredNorm() < Eigen::NumTraits<double>::dummy_precision());
  }

}
