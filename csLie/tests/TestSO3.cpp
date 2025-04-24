#include <cmath>
#include <csCore/math/Angles.hpp>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SO3_impl.hpp>
#include <csLie/Types.hpp>

namespace csLie {

/*
TEMPLATE_TEST_CASE("SO3 test - base", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  // create a identity rotation matrix
  SO3<S> const rI{SO3<S>::Identity()};

  // create a zero algebra
  SO3Alg<S> const d0{SO3Alg<S>::Zero()};

  // exp the zero algebra, it has to be the identity matrix
  SO3<S> const exp0{so3Exp(d0)};
  REQUIRE((exp0 - rI).isZero());

  // exp the 90deg (M_PI / 2), expect [0 -1; 1 0]
  SO3<S> const exp90{so3Exp<S>(csCore::math::degToRad(90.0))};
  SO3<S> exp90_gt{};
  exp90_gt << 0, -1, 1, 0;
  REQUIRE((exp90 - exp90_gt).isZero());

  S const PRECISION = Constants<S>::dummy_precision;

  // log of identity is 0
  REQUIRE_THAT(so3Log(rI), Catch::Matchers::WithinAbs(0.0, PRECISION));

  // log of exp90 is 90
  REQUIRE_THAT(so3Log(exp90), Catch::Matchers::WithinAbs(
                                  csCore::math::degToRad(90.0), PRECISION));

  // randomized tests
  Eigen::VectorX<S> randoms =
      Eigen::VectorX<S>::Random(100) * 4 * Constants<S>::pi;
  for (size_t i = 0; i < randoms.size(); i++) {
    REQUIRE(csCore::eigenUtils::isValidRotationMatrix2D(so3Exp(randoms(i))));
    // log(exp(d)) = d
    REQUIRE_THAT(
        csCore::math::angleDistance(so3Log(so3Exp(randoms(i))), randoms(i)),
        Catch::Matchers::WithinAbs(0.0, PRECISION));

    // exp(log(exp(d))) * exp(d)^T = I
    REQUIRE(
        (so3Exp(so3Log(so3Exp(randoms(i)))) * so3Exp(randoms(i)).transpose() -
         SO3<S>::Identity())
            .isZero());
  }
}

// TODO(sc) add a test with optimization in SO3

TEMPLATE_TEST_CASE("SO3 test - exp, delta wrong size", "[SO3]", float32_t,
                   float64_t) {
  using S = TestType;

  // so3Exp(Eigen::Matrix<S,3,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(so3Exp(Eigen::MatrixX<S>(1, 1)));
    REQUIRE_ASSERT(
        so3Exp(Eigen::MatrixX<S>(1, 2)));  // throws, wrong delta size
    REQUIRE_ASSERT(
        so3Exp(Eigen::MatrixX<S>(2, 1)));  // throws, wrong delta size

    REQUIRE_NO_ASSERT(so3Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 1, 1)));
    REQUIRE_NO_ASSERT(so3Exp(Eigen::VectorX<S>(10).segment(0, 1)));
    REQUIRE_ASSERT(so3Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(so3Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(so3Exp(Eigen::VectorX<S>(10).segment(0, 2)));
  }
}

TEMPLATE_TEST_CASE("SO3 test - log, input is not a rotation matrix", "[SO3]",
                   float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    // throws, wrong matrix size
    REQUIRE_ASSERT(so3Log(Eigen::MatrixX<S>(3, 3)));
    REQUIRE_ASSERT(so3Log(Eigen::MatrixX<S>(2, 3)));

    Eigen::MatrixX<S> m10x10 = Eigen::MatrixX<S>::Identity(10, 10);
    REQUIRE_NO_ASSERT(so3Log(
        m10x10.block(0, 0, 2, 2)));  // ok (is also a valid rotation matrix)!

    REQUIRE_ASSERT(so3Log(Eigen::MatrixX<S>(4, 4).block(0, 0, 3, 3)));
    REQUIRE_ASSERT(so3Log(Eigen::MatrixX<S>(4, 4).block(0, 0, 2, 1)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    REQUIRE_ASSERT(so3Log(
        Eigen::Matrix2<S>::Zero()));  // throws, zero is not a rotation matrix
  }
}
*/

TEMPLATE_TEST_CASE("SO3 test - from/to RPY", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  // 1 axis only
  {
    S a = csCore::math::degToRad(23.0);
    Eigen::Vector3<S> vz_gt(0, 0, a);
    Eigen::Vector3<S> vy_gt(0, a, 0);
    Eigen::Vector3<S> vx_gt(a, 0, 0);
    SO3<S> rz{rpyToSo3(vz_gt)};
    SO3<S> ry{rpyToSo3(vy_gt)};
    SO3<S> rx{rpyToSo3(vx_gt)};

    SO3<S> rz_gt;
    rz_gt = Eigen::AngleAxis<S>(a, Eigen::Vector3<S>::UnitZ());
    SO3<S> ry_gt;
    ry_gt = Eigen::AngleAxis<S>(a, Eigen::Vector3<S>::UnitY());
    SO3<S> rx_gt;
    rx_gt = Eigen::AngleAxis<S>(a, Eigen::Vector3<S>::UnitX());

    REQUIRE((rz_gt - rz).isZero());
    REQUIRE((ry_gt - ry).isZero());
    REQUIRE((rx_gt - rx).isZero());

    Eigen::Vector3<S> vz = so3ToRPY(rz);
    Eigen::Vector3<S> vy = so3ToRPY(ry);
    Eigen::Vector3<S> vx = so3ToRPY(rx);

    REQUIRE((vz_gt - vz).isZero());
    REQUIRE((vy_gt - vy).isZero());
    REQUIRE((vx_gt - vx).isZero());
  }
  // 3 axis
  {
    Eigen::Vector3<S> a(csCore::math::degToRad(23.0),
                        csCore::math::degToRad(57.0),
                        csCore::math::degToRad(82.0));
    SO3<S> r{rpyToSo3(a)};
    Eigen::Vector3<S> a_r{so3ToRPY(r)};

    SO3<S> r_e;
    r_e = Eigen::AngleAxis<S>(a(2), Eigen::Vector3<S>::UnitZ()) *
          Eigen::AngleAxis<S>(a(1), Eigen::Vector3<S>::UnitY()) *
          Eigen::AngleAxis<S>(a(0), Eigen::Vector3<S>::UnitX());

    REQUIRE((a - a_r).isZero());
    REQUIRE((r - r_e).isZero());
  }
}

// reference for RPY conversion
// https://eecs.qmul.ac.uk/~gslabaugh/publications/euler.pdf

TEMPLATE_TEST_CASE("SO3 test - Gymbal lock", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  // 90 deg
  {
    Eigen::Vector3<S> vy(csCore::math::degToRad(10.0),
                         csCore::math::degToRad(90.0),
                         csCore::math::degToRad(-20.0));
    SO3<S> r1{rpyToSo3(vy)};

    // gymbal on yaw - roll

    vy(0) = csCore::math::degToRad(30.0);
    vy(2) = csCore::math::degToRad(0.0);
    SO3<S> r2{rpyToSo3(vy)};

    vy(0) = csCore::math::degToRad(0.0);
    vy(2) = csCore::math::degToRad(-30.0);
    SO3<S> r3{rpyToSo3(vy)};

    vy(0) = csCore::math::degToRad(-10.0);
    vy(2) = csCore::math::degToRad(-40.0);
    SO3<S> r4{rpyToSo3(vy)};

    REQUIRE((r1 * r2.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
    REQUIRE((r1 * r3.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
    REQUIRE((r1 * r4.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
  }

  // -90 deg
  {
    Eigen::Vector3<S> vy(csCore::math::degToRad(10.0),
                         csCore::math::degToRad(-90.0),
                         csCore::math::degToRad(20.0));
    SO3<S> r1{rpyToSo3(vy)};

    // gymbal on yaw + roll

    vy(0) = csCore::math::degToRad(30.0);
    vy(2) = csCore::math::degToRad(0.0);
    SO3<S> r2{rpyToSo3(vy)};

    std::cout << "r1 " << std::endl << r1 << std::endl;
    std::cout << "r2 " << std::endl << r2 << std::endl;

    vy(0) = csCore::math::degToRad(0.0);
    vy(2) = csCore::math::degToRad(30.0);
    SO3<S> r3{rpyToSo3(vy)};

    vy(0) = csCore::math::degToRad(-10.0);
    vy(2) = csCore::math::degToRad(40.0);
    SO3<S> r4{rpyToSo3(vy)};

    REQUIRE((r1 * r2.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
    REQUIRE((r1 * r3.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
    REQUIRE((r1 * r4.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
  }
}

TEMPLATE_TEST_CASE("SO3 test - from RPY, input is not a 3d vector", "[SO3]",
                   float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    // throws, wrong matrix size
    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(3, 2)));
    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(1, 3)));
    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(2, 1)));

    Eigen::MatrixX<S> m10x10 = Eigen::MatrixX<S>::Identity(10, 10);
    REQUIRE_NO_ASSERT(rpyToSo3(
        m10x10.block(0, 0, 3, 1)));  // ok (is also a valid rotation matrix)!

    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(4, 4).block(0, 0, 3, 2)));
    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(4, 4).block(0, 0, 1, 3)));
  }
}

TEMPLATE_TEST_CASE("SO3 test - to RPY, input is not a rotation matrix", "[SO3]",
                   float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    // throws, wrong matrix size
    REQUIRE_ASSERT(so3ToRPY(Eigen::MatrixX<S>(3, 4)));
    REQUIRE_ASSERT(so3ToRPY(Eigen::MatrixX<S>(2, 3)));

    Eigen::MatrixX<S> m10x10 = Eigen::MatrixX<S>::Identity(10, 10);
    REQUIRE_NO_ASSERT(so3ToRPY(
        m10x10.block(0, 0, 3, 3)));  // ok (is also a valid rotation matrix)!

    REQUIRE_ASSERT(so3ToRPY(Eigen::MatrixX<S>(4, 4).block(0, 0, 2, 3)));
    REQUIRE_ASSERT(so3ToRPY(Eigen::MatrixX<S>(4, 4).block(0, 0, 2, 1)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    REQUIRE_ASSERT(so3ToRPY(
        Eigen::Matrix3<S>::Zero()));  // throws, zero is not a rotation matrix
  }
}

}  // namespace csLie
