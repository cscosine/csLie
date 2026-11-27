#include <cmath>
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
  REQUIRE(csCore::eigenUtils::isValidIsometryMatrix2D(exp0));

  REQUIRE((exp0.matrix() - rI.matrix()).isZero());

  // exp the 90deg (M_PI / 2) with a transl, expect [0 -1; 1 0] in rotation
  S const x = 10;
  S const y = 20;
  S const theta = csCore::math::degToRad(S(90.0));
  SE2<S> const exp90{se2Exp(SE2Alg<S>(x, y, theta))};
  REQUIRE(csCore::eigenUtils::isValidIsometryMatrix2D(exp90));

  SE2<S> exp90_gt{};

  S const sinTh_Th = S(1) / theta;
  S const m1cosTh_Th = S(1) / theta;
  S const tx = x * sinTh_Th - y * m1cosTh_Th;
  S const ty = x * m1cosTh_Th + y * sinTh_Th;

  exp90_gt.matrix() << 0, -1, tx, 1, 0, ty, 0, 0, 1;

  REQUIRE((exp90.matrix() - exp90_gt.matrix()).isZero());
}

TEMPLATE_TEST_CASE("SE2 test - from/to xyY", "[SE2]", float32_t, float64_t) {
  using S = TestType;

  S const x = 10;
  S const y = 20;
  S const theta = csCore::math::degToRad(S(90.0));

  Eigen::Vector3<S> xyY_gt(x, y, theta);

  SE2<S> transf_gt = SE2<S>::Identity();
  transf_gt.translation() << x, y;
  transf_gt.linear() << 0, -1, 1, 0;

  SE2<S> transf = xyYToSe2(xyY_gt);

  Eigen::Vector3<S> xyY = se2ToxyY(transf_gt);
  Eigen::Vector3<S> xyY_transf = se2ToxyY(transf_gt.matrix());

  REQUIRE((transf_gt.matrix() - transf.matrix()).isZero());

  REQUIRE((xyY_gt - xyY).isZero());
  REQUIRE((xyY_gt - xyY_transf).isZero());
}

TEMPLATE_TEST_CASE("SE2 test - log", "[SE2]", float32_t, float64_t) {
  using S = TestType;

  // log(exp(xi)) = xi for theta in (-pi, pi], exp(log(T)) = T
  for (S const theta : {S(0), S(1e-6), S(-0.3), S(1), S(-2), S(3)}) {
    SE2Alg<S> const xi(S(1.5), S(-0.25), theta);
    SE2<S> const T = se2Exp(xi);
    SE2Alg<S> const logT = se2Log(T);
    INFO("xi " << xi.transpose() << ", log(exp(xi)) " << logT.transpose());
    REQUIRE((logT - xi).isZero());
    REQUIRE((se2Exp(logT).matrix() - T.matrix()).isZero());
    REQUIRE((se2Log(T.matrix()) - logT).isZero());
  }

  // beyond pi the log is the equivalent angle in (-pi, pi], same rotation and translation
  SE2Alg<S> const xi(S(1), S(2), S(4));
  SE2Alg<S> const logT = se2Log(se2Exp(xi));
  REQUIRE(logT(2) <= Constants<S>::pi);
  REQUIRE((se2Exp(logT).matrix() - se2Exp(xi).matrix()).isZero());
}

// TODO(sc) add a test with optimization in SE2

TEMPLATE_TEST_CASE("SE2 test - exp, delta wrong size", "[SE2]", float32_t, float64_t) {
  using S = TestType;

  // se2Exp(Eigen::Matrix<S,4,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(se2Exp(Eigen::MatrixX<S>(3, 1)));
    REQUIRE_ASSERT(se2Exp(Eigen::MatrixX<S>(1, 2))); // throws, wrong size
    REQUIRE_ASSERT(se2Exp(Eigen::MatrixX<S>(2, 1))); // throws, wrong size

    REQUIRE_NO_ASSERT(se2Exp(Eigen::MatrixX<S>(4, 5).block(0, 0, 3, 1)));
    REQUIRE_NO_ASSERT(se2Exp(Eigen::VectorX<S>(10).segment(0, 3)));
    REQUIRE_ASSERT(se2Exp(Eigen::MatrixX<S>(6, 10).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(se2Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(se2Exp(Eigen::VectorX<S>(10).segment(0, 2)));
  }
}

TEMPLATE_TEST_CASE("SE2 test - log, transf wrong", "[SE2]", float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(se2Log(Eigen::MatrixX<S>::Identity(3, 3)));
    REQUIRE_ASSERT(se2Log(Eigen::MatrixX<S>(1, 2))); // throws, wrong size
    REQUIRE_ASSERT(se2Log(Eigen::MatrixX<S>(2, 1))); // throws, wrong size

    Eigen::MatrixX<S> m44 = Eigen::MatrixX<S>::Identity(4, 4);
    REQUIRE_NO_ASSERT(se2Log(m44.block(0, 0, 3, 3)));
    REQUIRE_ASSERT(se2Log(Eigen::MatrixX<S>(6, 10).block(0, 0, 1, 2)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    SE2<S> tmp = SE2<S>::Identity();
    tmp(2, 2) = 0.0;
    REQUIRE_ASSERT(se2Log(tmp));
  }
}

TEMPLATE_TEST_CASE("SE2 test - xyYToSe2, xyY wrong size", "[SE2]", float32_t, float64_t) {
  using S = TestType;

  // xyYToSe2(Eigen::Matrix<S,4,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(xyYToSe2(Eigen::MatrixX<S>(3, 1)));
    REQUIRE_ASSERT(xyYToSe2(Eigen::MatrixX<S>(1, 2))); // throws, wrong size
    REQUIRE_ASSERT(xyYToSe2(Eigen::MatrixX<S>(2, 1))); // throws, wrong size

    REQUIRE_NO_ASSERT(xyYToSe2(Eigen::MatrixX<S>(4, 5).block(0, 0, 3, 1)));
    REQUIRE_NO_ASSERT(xyYToSe2(Eigen::VectorX<S>(10).segment(0, 3)));
    REQUIRE_ASSERT(xyYToSe2(Eigen::MatrixX<S>(6, 10).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(xyYToSe2(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(xyYToSe2(Eigen::VectorX<S>(10).segment(0, 2)));
  }
}

TEMPLATE_TEST_CASE("SE2 test - se2ToxyY, se2 wrong", "[SE2]", float32_t, float64_t) {
  using S = TestType;

  // se2ToxyY(Eigen::Matrix<S,4,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(se2ToxyY(Eigen::MatrixX<S>::Identity(3, 3)));
    REQUIRE_ASSERT(se2ToxyY(Eigen::MatrixX<S>(1, 2))); // throws, wrong size
    REQUIRE_ASSERT(se2ToxyY(Eigen::MatrixX<S>(2, 1))); // throws, wrong size

    Eigen::MatrixX<S> m44 = Eigen::MatrixX<S>::Identity(4, 4);
    REQUIRE_NO_ASSERT(se2ToxyY(m44.block(0, 0, 3, 3)));
    REQUIRE_ASSERT(se2ToxyY(Eigen::MatrixX<S>(6, 10).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(se2ToxyY(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    SE2<S> tmp = SE2<S>::Identity();
    tmp(2, 2) = 0.0;
    REQUIRE_ASSERT(se2ToxyY(tmp));
  }
}

} // namespace csLie
