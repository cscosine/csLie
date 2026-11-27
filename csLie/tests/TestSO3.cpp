#include <cmath>
#include <csCore/math/Angles.hpp>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SO3_impl.hpp>
#include <csLie/Types.hpp>

namespace csLie {

TEMPLATE_TEST_CASE("SO3 test - base", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  // create a identity rotation matrix
  SO3<S> const rI{SO3<S>::Identity()};

  // create a zero algebra
  SO3Alg<S> const d0{SO3Alg<S>::Zero()};

  // exp the zero algebra, it has to be the identity matrix
  SO3<S> const exp0{so3Exp(d0)};
  REQUIRE((exp0 - rI).isZero());

  // exp the 90deg on z (M_PI / 2), expect [0 -1 0; 1 0 0; 0 0 1]
  SO3Alg<S> const rotZ90(0, 0, csCore::math::degToRad(S(90.0)));
  SO3<S> const exp90{so3Exp(rotZ90)};
  SO3<S> exp90_gt{};
  exp90_gt << 0, -1, 0, 1, 0, 0, 0, 0, 1;
  REQUIRE((exp90 - exp90_gt).isZero());

  // log of identity is 0
  REQUIRE(so3Log(rI).isZero());

  // log of exp90 is 90
  REQUIRE((so3Log(exp90) - rotZ90).isZero());

  // randomized tests
  for (size_t i = 0; i < 100; i++) {
    SO3Alg<S> const randV = SO3Alg<S>::Random();

    SO3<S> const R = so3Exp(randV);

    REQUIRE(csCore::eigenUtils::isValidRotationMatrix3D(R));

    // log(exp(d)) = d
    SO3Alg<S> const logV = so3Log(R);

    SO3<S> const RlogV = so3Exp(logV);
    REQUIRE((RlogV - R).isZero());

    // exp(log(exp(d))) * exp(d)^T = I
    REQUIRE((so3Exp(so3Log(so3Exp(randV))) * so3Exp(randV).transpose() - SO3<S>::Identity()).isZero());

    // exp(v) = eye + skewSymm(v)*J
    SO3<S> const RJ = Eigen::Matrix<S, 3, 3>::Identity() + skewSymm(randV) * so3J(randV);
    REQUIRE(((R * RJ.transpose()).matrix() - SO3<S>::Identity()).isZero());

    // check J and Jinv
    SO3J<S> const J{so3J(randV)};
    SO3Jinv<S> const Jinv{so3Jinv(randV)};
    INFO("J\n" << J << "\n\nJinv\n" << Jinv << "\n\nexpected eye\n" << J * Jinv);
    REQUIRE((J * Jinv - Eigen::Matrix<S, 3, 3>::Identity()).isZero());
  }

  // TODO(sc) cover log and exp limit cases
}

// TODO(sc) add a test with optimization in SO3
// TODO(sc) add a test with optimization in SO3 using J and Jinv

TEMPLATE_TEST_CASE("SO3 test - from/to RPY", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  // 1 axis only
  {
    S a = csCore::math::degToRad(S(23.0));
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

    SO3<S> rvz{rpyToSo3(vz)};
    SO3<S> rvy{rpyToSo3(vy)};
    SO3<S> rvx{rpyToSo3(vx)};

    REQUIRE((rz_gt - rvz).isZero());
    REQUIRE((ry_gt - rvy).isZero());
    REQUIRE((rx_gt - rvx).isZero());
  }
  // 3 axis
  {
    Eigen::Vector3<S> a(csCore::math::degToRad(S(23.0)), csCore::math::degToRad(S(57.0)),
                        csCore::math::degToRad(S(82.0)));
    SO3<S> r{rpyToSo3(a)};
    Eigen::Vector3<S> a_r{so3ToRPY(r)};
    SO3<S> r_a_r{rpyToSo3(a_r)};

    SO3<S> r_e;
    r_e = Eigen::AngleAxis<S>(a(2), Eigen::Vector3<S>::UnitZ()) *
          Eigen::AngleAxis<S>(a(1), Eigen::Vector3<S>::UnitY()) * Eigen::AngleAxis<S>(a(0), Eigen::Vector3<S>::UnitX());

    REQUIRE((r_a_r - r).isZero());
    REQUIRE((r - r_e).isZero());
  }
}

TEMPLATE_TEST_CASE("SO3 test - Gymbal lock", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  // 90 deg
  {
    Eigen::Vector3<S> vy(csCore::math::degToRad(S(10.0)), csCore::math::degToRad(S(90.0)),
                         csCore::math::degToRad(S(-20.0)));
    SO3<S> r1{rpyToSo3(vy)};

    // gymbal on yaw - roll

    vy(0) = csCore::math::degToRad(S(30.0));
    vy(2) = csCore::math::degToRad(S(0.0));
    SO3<S> r2{rpyToSo3(vy)};

    vy(0) = csCore::math::degToRad(S(0.0));
    vy(2) = csCore::math::degToRad(S(-30.0));
    SO3<S> r3{rpyToSo3(vy)};

    vy(0) = csCore::math::degToRad(S(-10.0));
    vy(2) = csCore::math::degToRad(S(-40.0));
    SO3<S> r4{rpyToSo3(vy)};

    REQUIRE((r1 * r2.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
    REQUIRE((r1 * r3.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
    REQUIRE((r1 * r4.transpose() - Eigen::Matrix3<S>::Identity()).isZero());

    Eigen::Vector3<S> v1{so3ToRPY(r1)};
    SO3<S> rv1{rpyToSo3(v1)};
    REQUIRE((rv1 * r1.transpose() - Eigen::Matrix3<S>::Identity()).isZero());

    Eigen::Vector3<S> v2{so3ToRPY(r2)};
    SO3<S> rv2{rpyToSo3(v2)};
    REQUIRE((rv2 * r2.transpose() - Eigen::Matrix3<S>::Identity()).isZero());

    Eigen::Vector3<S> v3{so3ToRPY(r3)};
    SO3<S> rv3{rpyToSo3(v3)};
    REQUIRE((rv3 * r3.transpose() - Eigen::Matrix3<S>::Identity()).isZero());

    Eigen::Vector3<S> v4{so3ToRPY(r4)};
    SO3<S> rv4{rpyToSo3(v4)};
    REQUIRE((rv4 * r4.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
  }

  // -90 deg
  {
    Eigen::Vector3<S> vy(csCore::math::degToRad(S(10.0)), csCore::math::degToRad(S(-90.0)),
                         csCore::math::degToRad(S(20.0)));
    SO3<S> r1{rpyToSo3(vy)};

    // gymbal on yaw + roll

    vy(0) = csCore::math::degToRad(S(30.0));
    vy(2) = csCore::math::degToRad(S(0.0));
    SO3<S> r2{rpyToSo3(vy)};

    vy(0) = csCore::math::degToRad(S(0.0));
    vy(2) = csCore::math::degToRad(S(30.0));
    SO3<S> r3{rpyToSo3(vy)};

    vy(0) = csCore::math::degToRad(S(-10.0));
    vy(2) = csCore::math::degToRad(S(40.0));
    SO3<S> r4{rpyToSo3(vy)};

    REQUIRE((r1 * r2.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
    REQUIRE((r1 * r3.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
    REQUIRE((r1 * r4.transpose() - Eigen::Matrix3<S>::Identity()).isZero());

    Eigen::Vector3<S> v1{so3ToRPY(r1)};
    SO3<S> rv1{rpyToSo3(v1)};
    REQUIRE((rv1 * r1.transpose() - Eigen::Matrix3<S>::Identity()).isZero());

    Eigen::Vector3<S> v2{so3ToRPY(r2)};
    SO3<S> rv2{rpyToSo3(v2)};
    REQUIRE((rv2 * r2.transpose() - Eigen::Matrix3<S>::Identity()).isZero());

    Eigen::Vector3<S> v3{so3ToRPY(r3)};
    SO3<S> rv3{rpyToSo3(v3)};
    REQUIRE((rv3 * r3.transpose() - Eigen::Matrix3<S>::Identity()).isZero());

    Eigen::Vector3<S> v4{so3ToRPY(r4)};
    SO3<S> rv4{rpyToSo3(v4)};
    REQUIRE((rv4 * r4.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
  }
}

TEMPLATE_TEST_CASE("SO3 test - SkewSymm", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  Eigen::Matrix3<S> s{skewSymm(Eigen::Vector3<S>(S(2.4), S(1.2), S(3.2)))};
  REQUIRE(s.diagonal().isZero(0));
  REQUIRE((s + s.transpose()).isZero(0));
}

TEMPLATE_TEST_CASE("SO3 test - exp, delta wrong size", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  // so3Exp(Eigen::Matrix<S,4,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(so3Exp(Eigen::MatrixX<S>(3, 1)));
    REQUIRE_ASSERT(so3Exp(Eigen::MatrixX<S>(1, 2))); // throws, wrong delta size
    REQUIRE_ASSERT(so3Exp(Eigen::MatrixX<S>(2, 1))); // throws, wrong delta size

    REQUIRE_NO_ASSERT(so3Exp(Eigen::MatrixX<S>(5, 7).block(0, 0, 3, 1)));
    REQUIRE_NO_ASSERT(so3Exp(Eigen::VectorX<S>(10).segment(0, 3)));
    REQUIRE_ASSERT(so3Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(so3Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(so3Exp(Eigen::VectorX<S>(10).segment(0, 2)));
  }
}

TEMPLATE_TEST_CASE("SO3 test - SO3J, delta wrong size", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  // so3J(Eigen::Matrix<S,4,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(so3J(Eigen::MatrixX<S>(3, 1)));
    REQUIRE_ASSERT(so3J(Eigen::MatrixX<S>(1, 2))); // throws, wrong delta size
    REQUIRE_ASSERT(so3J(Eigen::MatrixX<S>(2, 1))); // throws, wrong delta size

    REQUIRE_NO_ASSERT(so3J(Eigen::MatrixX<S>(5, 7).block(0, 0, 3, 1)));
    REQUIRE_NO_ASSERT(so3J(Eigen::VectorX<S>(10).segment(0, 3)));
    REQUIRE_ASSERT(so3J(Eigen::MatrixX<S>(2, 2).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(so3J(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(so3J(Eigen::VectorX<S>(10).segment(0, 2)));
  }
}

TEMPLATE_TEST_CASE("SO3 test - so3Jinv, delta wrong size", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  // so3J(Eigen::Matrix<S,4,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(so3Jinv(Eigen::MatrixX<S>(3, 1)));
    REQUIRE_ASSERT(so3Jinv(Eigen::MatrixX<S>(1, 2))); // throws, wrong delta size
    REQUIRE_ASSERT(so3Jinv(Eigen::MatrixX<S>(2, 1))); // throws, wrong delta size

    REQUIRE_NO_ASSERT(so3Jinv(Eigen::MatrixX<S>(5, 7).block(0, 0, 3, 1)));
    REQUIRE_NO_ASSERT(so3Jinv(Eigen::VectorX<S>(10).segment(0, 3)));
    REQUIRE_ASSERT(so3Jinv(Eigen::MatrixX<S>(2, 2).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(so3Jinv(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(so3Jinv(Eigen::VectorX<S>(10).segment(0, 2)));
  }
}

TEMPLATE_TEST_CASE("SO3 test - log, input is not a rotation matrix", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    // throws, wrong matrix size
    REQUIRE_ASSERT(so3Log(Eigen::MatrixX<S>(2, 2)));
    REQUIRE_ASSERT(so3Log(Eigen::MatrixX<S>(4, 5)));

    Eigen::MatrixX<S> m10x10 = Eigen::MatrixX<S>::Identity(10, 10);
    REQUIRE_NO_ASSERT(so3Log(m10x10.block(0, 0, 3, 3))); // ok (is also a valid rotation matrix)!

    REQUIRE_ASSERT(so3Log(Eigen::MatrixX<S>(4, 4).block(0, 0, 4, 3)));
    REQUIRE_ASSERT(so3Log(Eigen::MatrixX<S>(4, 4).block(0, 0, 2, 1)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    REQUIRE_ASSERT(so3Log(Eigen::Matrix3<S>::Zero())); // throws, zero is not a rotation matrix
  }
}

TEMPLATE_TEST_CASE("SO3 test - from RPY, input is not a 3d vector", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    // throws, wrong matrix size
    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(3, 2)));
    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(1, 3)));
    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(2, 1)));

    Eigen::MatrixX<S> m10x10 = Eigen::MatrixX<S>::Identity(10, 10);
    REQUIRE_NO_ASSERT(rpyToSo3(m10x10.block(0, 0, 3, 1)));

    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(4, 4).block(0, 0, 3, 2)));
    REQUIRE_ASSERT(rpyToSo3(Eigen::MatrixX<S>(4, 4).block(0, 0, 1, 3)));
  }
}

TEMPLATE_TEST_CASE("SO3 test - to RPY, input is not a rotation matrix", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    // throws, wrong matrix size
    REQUIRE_ASSERT(so3ToRPY(Eigen::MatrixX<S>(3, 4)));
    REQUIRE_ASSERT(so3ToRPY(Eigen::MatrixX<S>(2, 3)));

    Eigen::MatrixX<S> m10x10 = Eigen::MatrixX<S>::Identity(10, 10);
    REQUIRE_NO_ASSERT(so3ToRPY(m10x10.block(0, 0, 3, 3))); // ok (is also a valid rotation matrix)!

    REQUIRE_ASSERT(so3ToRPY(Eigen::MatrixX<S>(4, 4).block(0, 0, 2, 3)));
    REQUIRE_ASSERT(so3ToRPY(Eigen::MatrixX<S>(4, 4).block(0, 0, 2, 1)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    REQUIRE_ASSERT(so3ToRPY(Eigen::Matrix3<S>::Zero())); // throws, zero is not a rotation matrix
  }
}

TEMPLATE_TEST_CASE("SO3 test - skewSymm, input is not a 3d vector", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    // throws, wrong matrix size
    REQUIRE_ASSERT(skewSymm(Eigen::MatrixX<S>(3, 2)));
    REQUIRE_ASSERT(skewSymm(Eigen::MatrixX<S>(1, 3)));
    REQUIRE_ASSERT(skewSymm(Eigen::MatrixX<S>(2, 1)));

    Eigen::MatrixX<S> m10x10 = Eigen::MatrixX<S>::Identity(10, 10);
    REQUIRE_NO_ASSERT(skewSymm(m10x10.block(0, 0, 3, 1)));

    REQUIRE_ASSERT(skewSymm(Eigen::MatrixX<S>(4, 4).block(0, 0, 3, 2)));
    REQUIRE_ASSERT(skewSymm(Eigen::MatrixX<S>(4, 4).block(0, 0, 1, 3)));
  }
}

} // namespace csLie
