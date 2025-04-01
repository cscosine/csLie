#include <csCore/testUtils/AssertCatch2.hpp>

#include <csLie/Types.hpp>
#include <csLie/SO2_impl.hpp>
#include <csLie/Constants.hpp>

#include <csLie/tests/TestUtils.hpp>

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
    SO2<S> const exp90{ so2Exp<S>(Constants<S>::pi / static_cast<S>(2.0)) };
    SO2<S> exp90_gt{};
    exp90_gt << 
      0, -1,
      1, 0;
    REQUIRE((exp90 - exp90_gt).isZero());

  }

  TEMPLATE_TEST_CASE("SO2 test - delta wrong size", "[SO2]", float32_t, float64_t) 
  {
    using S = TestType;

    // so2Exp(Eigen::Matrix<S,3,1>()); // does not compile
    if(csLieAsserts::isAssertDebugEnabled()) {
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