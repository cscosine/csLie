#include <cmath>
#include <csCore/eigenUtils/RotationMatrix_impl.hpp>
#include <csCore/math/Angles.hpp>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SE2_impl.hpp>
#include <csLie/Types.hpp>

namespace csLie {

TEMPLATE_TEST_CASE("SE2 test - base", "[SE2]", float32_t, float64_t) {
  using S = TestType;

  // create a identity matrix
  SE2<S> const rI{SE2<S>::Identity()};

  // create a zero algebra
  SE2Alg<S> const d0{SE2Alg<S>::Zero()};

  // exp the zero algebra, it has to be the identity matrix
  SE2<S> const exp0{se2Exp(d0)};
  REQUIRE((exp0.matrix() - rI.matrix()).isZero());

  // exp the 90deg (M_PI / 2) with zero transl, expect [0 -1; 1 0]
  SE2<S> const exp90{se2Exp(SE2Alg<S>(0, 0, csCore::math::degToRad(90.0)))};
  SE2<S> exp90_gt{};
  exp90_gt.linear() << 0, -1, 1, 0;
  REQUIRE((exp90.matrix() - exp90_gt.matrix()).isZero());
}

// TODO(sc) add test of log in SE2, after adding log
// TODO(sc) add a test with optimization in SE2

TEMPLATE_TEST_CASE("SE2 test - exp, delta wrong size", "[SE2]", float32_t,
                   float64_t) {
  using S = TestType;

  // se2Exp(Eigen::Matrix<S,3,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(se2Exp(Eigen::MatrixX<S>(3, 1)));
    REQUIRE_ASSERT(
        se2Exp(Eigen::MatrixX<S>(1, 2)));  // throws, wrong delta size
    REQUIRE_ASSERT(
        se2Exp(Eigen::MatrixX<S>(2, 1)));  // throws, wrong delta size

    REQUIRE_NO_ASSERT(se2Exp(Eigen::MatrixX<S>(4, 5).block(0, 0, 3, 1)));
    REQUIRE_NO_ASSERT(se2Exp(Eigen::VectorX<S>(10).segment(0, 3)));
    REQUIRE_ASSERT(se2Exp(Eigen::MatrixX<S>(6, 10).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(se2Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(se2Exp(Eigen::VectorX<S>(10).segment(0, 2)));
  }
}

}  // namespace csLie
