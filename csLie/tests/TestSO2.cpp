#include <cmath>
#include <csCore/testUtils/AssertCatch2.hpp>

#include <csCore/eigenUtils/RotationMatrix_impl.hpp>

#include <csLie/Types.hpp>
#include <csLie/SO2_impl.hpp>
#include <csLie/Constants.hpp>

#include <csCore/math/Angles.hpp>

namespace csLie {
  

  TEMPLATE_TEST_CASE("SO2 test - base", "[SO2]", float32_t, float64_t) {
    using S = TestType;

    // create a identity rotation matrix
    SO2<S> const rI{ SO2<S>::Identity()};

    // create a zero algebra
    SO2Alg<S> const d0{ SO2Alg<S>::Zero() };

    // exp the zero algebra, it has to be the identity matrix
    SO2<S> const exp0{ so2Exp(d0) };
    REQUIRE((exp0 - rI).isZero());

    // exp the 90deg (M_PI / 2), expect [0 -1; 1 0]
    SO2<S> const exp90{ so2Exp<S>(csCore::math::degToRad(90.0)) };
    SO2<S> exp90_gt{};
    exp90_gt << 
      0, -1,
      1, 0;
    REQUIRE((exp90 - exp90_gt).isZero());

    S const PRECISION = Constants<S>::dummy_precision;

    // log of identity is 0
    REQUIRE_THAT(so2Log(rI), Catch::Matchers::WithinAbs(0.0, PRECISION));

    // log of exp90 is 90
    REQUIRE_THAT(so2Log(exp90), Catch::Matchers::WithinAbs(csCore::math::degToRad(90.0), PRECISION));


    // randomized tests
    Eigen::VectorX<S> randoms = Eigen::VectorX<S>::Random(100) * 4 * Constants<S>::pi;
    for(size_t i = 0; i < randoms.size();i++) {
      REQUIRE(csCore::eigenUtils::isValidRotationMatrix2D(so2Exp(randoms(i))));
      // log(exp(d)) = d
      REQUIRE_THAT(csCore::math::angleDistance(so2Log(so2Exp(randoms(i))),randoms(i)), Catch::Matchers::WithinAbs(0.0, PRECISION));

      // exp(log(exp(d))) * exp(d)^T = I
      REQUIRE((so2Exp(so2Log(so2Exp(randoms(i)))) * so2Exp(randoms(i)).transpose() - SO2<S>::Identity()).isZero());
    }

  }

  TEMPLATE_TEST_CASE("SO2 test - exp, delta wrong size", "[SO2]", float32_t, float64_t) 
  {
    using S = TestType;

    // so2Exp(Eigen::Matrix<S,3,1>()); // does not compile
    if(csLie::Asserts::isAssertStandardEnabled()) {
      REQUIRE_NO_ASSERT(so2Exp(Eigen::MatrixX<S>(1,1))); 
      REQUIRE_ASSERT(so2Exp(Eigen::MatrixX<S>(1,2))); // throws, wrong delta size
      REQUIRE_ASSERT(so2Exp(Eigen::MatrixX<S>(2,1))); // throws, wrong delta size

      REQUIRE_NO_ASSERT(so2Exp(Eigen::MatrixX<S>(2,2).block(0,0,1,1))); 
      REQUIRE_NO_ASSERT(so2Exp(Eigen::VectorX<S>(10).segment(0,1))); 
      REQUIRE_ASSERT(so2Exp(Eigen::MatrixX<S>(2,2).block(0,0,1,2))); 
      REQUIRE_ASSERT(so2Exp(Eigen::MatrixX<S>(2,2).block(0,0,2,1))); 
      REQUIRE_ASSERT(so2Exp(Eigen::VectorX<S>(10).segment(0,2))); 
    }
    
  }

}