#include <cmath>
#include <csCore/math/Angles.hpp>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SE3_impl.hpp>
#include <csLie/Types.hpp>

namespace csLie {

TEMPLATE_TEST_CASE("SE3 test - base", "[SE3]", float32_t, float64_t) {
  using S = TestType;

  // create a identity transformation matrix
  SE3<S> const rI{SE3<S>::Identity()};

  // create a zero algebra
  SE3Alg<S> const d0{SE3Alg<S>::Zero()};

  // exp the zero algebra, it has to be the identity matrix
  SE3<S> const exp0{se3Exp(d0)};
  REQUIRE((exp0.matrix() - rI.matrix()).isZero());

  S const PRECISION = Constants<S>::dummy_precision;

  // log of identity is 0
  REQUIRE(se3Log(rI).isZero());

  // randomized tests
  for (size_t i = 0; i < 100; i++) {
    SE3Alg<S> const randV = SE3Alg<S>::Random() * 4 * Constants<S>::pi;

    SE3<S> const T = se3Exp(randV);
    REQUIRE(csCore::eigenUtils::isValidIsometryMatrix3D(T));

    // log(exp(d)) = d
    SE3Alg<S> const logV = se3Log(T);

    SE3<S> const TlogV = se3Exp(logV);
    REQUIRE(csCore::eigenUtils::isValidIsometryMatrix3D(TlogV));
    REQUIRE((TlogV.matrix() - T.matrix()).isZero());

    // exp(log(exp(d))) * exp(d)^T = I
    Eigen::Transform<S, 3, Eigen::Isometry> check{se3Exp(randV)};
    REQUIRE(((se3Exp(se3Log(se3Exp(randV))) * check.inverse()).matrix() -
             SE3<S>::Identity().matrix())
                .isZero());
  }

  // TODO(sc) cover log and exp limit cases
}

TEMPLATE_TEST_CASE("SE3 test - adjoint", "[SE3]", float32_t, float64_t) {
  using S = TestType;
  for (size_t i = 0; i < 100; i++) {
    SE3Alg<S> const v = SE3Alg<S>::Random();
    SE3Alg<S> const t = SE3Alg<S>::Random();

    SE3<S> const V = se3Exp(v);
    REQUIRE(csCore::eigenUtils::isValidIsometryMatrix3D(V));
    SE3<S> const T = se3Exp(t);
    REQUIRE(csCore::eigenUtils::isValidIsometryMatrix3D(T));

    SE3Adj<S> const AdT = se3Adjoint(T);

    SE3<S> const TV = T * V;
    SE3<S> const expAdTVT = se3Exp(AdT * v) * T;

    SE3<S> expected_eye =
        expAdTVT * Eigen::Transform<S, 3, Eigen::Isometry>(TV).inverse();
    REQUIRE((expected_eye.matrix() - SE3<S>::Identity().matrix()).isZero());
  }
}

// TODO(sc) add a test with optimization in SE3

TEMPLATE_TEST_CASE("SE3 test - exp, delta wrong size", "[SE3]", float32_t,
                   float64_t) {
  using S = TestType;

  // se3Exp(Eigen::Matrix<S,3,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(se3Exp(Eigen::MatrixX<S>(6, 1)));
    REQUIRE_ASSERT(
        se3Exp(Eigen::MatrixX<S>(1, 2)));  // throws, wrong delta size
    REQUIRE_ASSERT(
        se3Exp(Eigen::MatrixX<S>(2, 1)));  // throws, wrong delta size

    REQUIRE_NO_ASSERT(se3Exp(Eigen::MatrixX<S>(9, 7).block(0, 0, 6, 1)));
    REQUIRE_NO_ASSERT(se3Exp(Eigen::VectorX<S>(10).segment(0, 6)));
    REQUIRE_ASSERT(se3Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(se3Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(se3Exp(Eigen::VectorX<S>(10).segment(0, 2)));
  }
}

TEMPLATE_TEST_CASE("SE3 test - log, input is not a transformation matrix",
                   "[SE3]", float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    // throws, wrong matrix size
    REQUIRE_ASSERT(se3Log(Eigen::MatrixX<S>(2, 2)));
    REQUIRE_ASSERT(se3Log(Eigen::MatrixX<S>(4, 5)));

    Eigen::MatrixX<S> m10x10 = Eigen::MatrixX<S>::Identity(10, 10);
    REQUIRE_NO_ASSERT(se3Log(m10x10.block(
        0, 0, 4, 4)));  // ok (is also a valid transformation matrix)!

    REQUIRE_ASSERT(se3Log(Eigen::MatrixX<S>(4, 4).block(0, 0, 4, 3)));
    REQUIRE_ASSERT(se3Log(Eigen::MatrixX<S>(4, 4).block(0, 0, 2, 1)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    REQUIRE_ASSERT(se3Log(Eigen::Matrix4<S>::Zero()));  // throws, zero is not a
                                                        // transformation matrix
  }
}

TEMPLATE_TEST_CASE("SE3 test - adjoint, input is not a transformation matrix",
                   "[SE3]", float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    // throws, wrong matrix size
    REQUIRE_ASSERT(se3Adjoint(Eigen::MatrixX<S>(2, 2)));
    REQUIRE_ASSERT(se3Adjoint(Eigen::MatrixX<S>(4, 5)));

    Eigen::MatrixX<S> m10x10 = Eigen::MatrixX<S>::Identity(10, 10);
    REQUIRE_NO_ASSERT(se3Adjoint(m10x10.block(
        0, 0, 4, 4)));  // ok (is also a valid transformation matrix)!

    REQUIRE_ASSERT(se3Adjoint(Eigen::MatrixX<S>(4, 4).block(0, 0, 4, 3)));
    REQUIRE_ASSERT(se3Adjoint(Eigen::MatrixX<S>(4, 4).block(0, 0, 2, 1)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    REQUIRE_ASSERT(
        se3Adjoint(Eigen::Matrix4<S>::Zero()));  // throws, zero is not a
                                                 // transformation matrix
  }
}

}  // namespace csLie
