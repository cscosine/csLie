#include <Eigen/Dense> // for determinant
#include <catch2/catch_all.hpp>
#include <csCore/eigenUtils/RotationMatrix_impl.hpp>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/S2_impl.hpp>
#include <csLie/Types.hpp>

namespace csLie {

TEMPLATE_TEST_CASE("S2 test - verify assumption atan2(0,0) = 0", "[S2]", float32_t, float64_t) {
  // in s2toZAlignedSO3 we assume that atan2(0,0) = 0. This might be platform
  // dependent. In case this test fail, adapt the logic in s2toZAlignedSO3 to
  // handle this edge case
  REQUIRE(std::atan2<TestType>(0, 0) == 0.0);
}

TEMPLATE_TEST_CASE("S2 test - base", "[S2]", float32_t, float64_t) {
  using S = TestType;
  S const PRECISION = Constants<S>::dummy_precision;

  // create a unit vector v
  S2<S> const v{S2<S>::Random().normalized()};

  // verify that v is a unit vector (test self-sanity)
  REQUIRE_THAT(v.norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));

  // create rotation matrix from v
  SO3<S> const Rv = s2toZAlignedSO3(v);
  // check is a valid 3d rot matrix
  REQUIRE(csCore::eigenUtils::isValidRotationMatrix3D(Rv));
  // and its first column is v
  REQUIRE((Rv.col(0) - v).isZero());

  // create a random s2 algebra vector delta
  S2Alg<S> const delta{S2Alg<S>::Random()};

  // create a unit vector vd as exp of delta
  S2<S> const vd = s2Exp(delta);
  // TODO(sc) check edge cases (taylor expansion for small numbers)
  REQUIRE_THAT(vd.norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));

  // compute oplus(v,delta) explicitely in voe
  S2<S> const voe{Rv * vd};
  // verify voe is a unit vector
  REQUIRE_THAT(voe.norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));

  // compute oplus(v,delta) with s2oplus in vo
  S2<S> const vo{s2oplus(v, delta)};
  // verify vo is a unit vector
  REQUIRE_THAT(vo.norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));
  // verify vo is equivalent to voe
  REQUIRE((vo - voe).isZero());

  // compute the jacobian of oplus(v,delta) wrt delta, evalauted in delta = 0
  s2oplusJac<S> const J = s2_d_oplus_d_delta_0(v);
  // verify is equivalent to the first two columns of Rv
  REQUIRE((J - Rv.template rightCols<2>()).isZero());
}

TEMPLATE_TEST_CASE("S2 test - plane fitting", "[S2]", float32_t, float64_t) {
  // repeat the test in 4 different modes
  bool const pointsNoisedFlag = GENERATE(false, true);
  bool const perturbedInitialization = GENERATE(false, true);
  CAPTURE(pointsNoisedFlag, perturbedInitialization);

  using S = TestType;
  S const PRECISION = Constants<S>::dummy_precision;

  // generate a random rotation matrix
  Eigen::Transform<S, 3, Eigen::Isometry> Rt{Eigen::Transform<S, 3, Eigen::Isometry>::Identity()};
  Rt.translation().setRandom();
  Rt.linear() =
      Eigen::Matrix<S, 3, 3>(Eigen::AngleAxis<S>(Eigen::Matrix<S, 1, 1>::Random()(0), Eigen::Vector3<S>::UnitX()) *
                             Eigen::AngleAxis<S>(Eigen::Matrix<S, 1, 1>::Random()(0), Eigen::Vector3<S>::UnitY()) *
                             Eigen::AngleAxis<S>(Eigen::Matrix<S, 1, 1>::Random()(0), Eigen::Vector3<S>::UnitZ()));

  // generate a random set of points on the Z=0 plane defined by Rt
  int const N = 40;
  Eigen::Matrix<S, 3, Eigen::Dynamic> points(3, N);

  Eigen::Vector3<S> const pointCenter = Rt.translation();
  Eigen::Vector3<S> const normalGT = Rt.linear().col(2);

  for (int i = 0; i < N; i++) {
    // generate a random point with Z=0
    Eigen::Vector3<S> pr = Eigen::Vector3<S>::Random();
    pr.z() = 0;
    // and transform in the world frame
    points.col(i) = Rt * pr;
    // check that the point is effectively on the plane as expected
    // i.e., a vector from the center of the reference system to the point is
    // orthogonal to the plane normal
    S const ep = (pointCenter - points.col(i)).dot(normalGT);
    REQUIRE_THAT(ep, Catch::Matchers::WithinAbs(0.0, PRECISION));
  }

  Eigen::Matrix<S, 3, Eigen::Dynamic> pointsNoised = points;
  if (pointsNoisedFlag) {
    // eventually, add noise to points
    pointsNoised += Eigen::Matrix<S, 3, Eigen::Dynamic>::Random(3, N) * 0.02;
  }

  Eigen::Vector3<S> estNormal = normalGT;
  if (perturbedInitialization) {
    // eventually, add noise to initial estimation
    estNormal = s2oplus(estNormal, Eigen::Vector2<S>::Random() * 0.01);
    estNormal.normalize();
  }
  Eigen::Vector3<S> const estNormalInitial = estNormal;

  // note: the plane in homogeneous coordinate is
  // pi = [n d]
  // where
  // - n is the normal
  // - d = -c^t * n
  // where c is the center of the reference system
  // proof: pi dot c = 0, allow to compute d
  // the point plane error e is  e = p dot pi = [p 1] dot [n d] = p^t n + d =
  // p^t n -c^t n = (p - c)^t * n the derivative d_e_d_delta evaluated in zero,
  // computed on (p - c)^t * oplus(n, delta), is (p - c)^t *
  // s2_d_oplus_d_delta_0(n)

  // perform a standard Gauss Newton optimization loop
  S chi2 = 0;
  S initial_chi2 = 0;
  for (int it = 0; it < 10; it++) {
    // initialize pseudo hessian H matrix and b vector to zero
    Eigen::Matrix2<S> H = Eigen::Matrix2<S>::Zero();
    Eigen::Vector2<S> b = Eigen::Vector2<S>::Zero();
    chi2 = 0;
    // per each point
    for (int i = 0; i < N; i++) {
      // compute the point plane error
      Eigen::Vector3<S> const pd = (pointsNoised.col(i) - pointCenter);
      S const e = pd.dot(estNormal);
      // and the jacobian
      Eigen::Matrix<S, 1, 2> d_e_d_delta = pd.transpose() * s2_d_oplus_d_delta_0(estNormal);

      // cumulate terms
      S const c2 = e * e;
      H += d_e_d_delta.transpose() * d_e_d_delta;
      b += d_e_d_delta.transpose() * e;
      chi2 += c2;
    }

    // store initial chi2 for testing purposes
    if (it == 0) {
      initial_chi2 = chi2;
    }

    // compute the solution
    Eigen::LDLT<Eigen::Matrix2<S>> ldlt;
    ldlt.compute(H);
    Eigen::Vector2<S> x = ldlt.solve(-b);
    // update the estimate
    estNormal = s2oplus(estNormal, x);
    estNormal.normalize();
  }

  // final checks
  if (!pointsNoisedFlag) {
    // if no noise, expect to converge exactly at zero
    REQUIRE_THAT(chi2, Catch::Matchers::WithinAbs(0.0, PRECISION));
    REQUIRE((estNormal - normalGT).isZero());
    // if no noise, estimate of normal is exact, even when the start is
    // perturbed
    S const dotFinal = estNormal.dot(normalGT);
    REQUIRE_THAT(dotFinal, Catch::Matchers::WithinAbs(1.0, PRECISION));
  } else {
    // if noise we can just ask that the chi2 is decreasing
    REQUIRE(chi2 < initial_chi2);
  }

  if (perturbedInitialization && !pointsNoisedFlag) {
    // if the init is perturbed, but no noise on points, we expect that the
    // final estimation is closer to GT
    S const dotInitial = estNormalInitial.dot(normalGT);
    S const dotFinal = estNormal.dot(normalGT);
    REQUIRE(dotFinal > dotInitial);
  }
}

TEMPLATE_TEST_CASE("S2 test - assertions", "[S2]", float32_t, float64_t) {
  using S = TestType;

  // s2Exp(Eigen::Matrix<S,3,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(s2Exp(Eigen::MatrixX<S>(2, 1)));
    REQUIRE_ASSERT(s2Exp(Eigen::MatrixX<S>(1, 2))); // throws, wrong delta size
    REQUIRE_ASSERT(s2Exp(Eigen::MatrixX<S>(2, 2))); // throws, wrong delta size
    REQUIRE_ASSERT(s2Exp(Eigen::MatrixX<S>(3, 2))); // throws, wrong delta size

    REQUIRE_NO_ASSERT(s2Exp(Eigen::MatrixX<S>(2, 2).block(0, 0, 2, 1)));
    REQUIRE_NO_ASSERT(s2Exp(Eigen::VectorX<S>(10).segment(0, 2)));
    REQUIRE_ASSERT(s2Exp(Eigen::MatrixX<S>(4, 3).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(s2Exp(Eigen::MatrixX<S>(4, 3).block(0, 0, 3, 1)));
    REQUIRE_ASSERT(s2Exp(Eigen::MatrixX<S>(4, 3).block(0, 0, 3, 3)));
    REQUIRE_ASSERT(s2Exp(Eigen::VectorX<S>(10).segment(0, 3)));
    REQUIRE_ASSERT(s2Exp(Eigen::VectorX<S>(10).segment(0, 1)));
  }

  // s2toZAlignedSO3(Eigen::Matrix<S,2,1>()); // does not compile
  // s2_d_oplus_d_delta_0(Eigen::Matrix<S,2,1>()); // does not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    Eigen::MatrixX<S> m31(3, 1);
    m31.setZero();
    m31(0, 0) = 1;
    REQUIRE_NO_ASSERT(s2toZAlignedSO3(m31));
    REQUIRE_NO_ASSERT(s2_d_oplus_d_delta_0(m31));

    Eigen::VectorX<S> v3(3);
    v3.setZero();
    v3(2) = 1;
    REQUIRE_NO_ASSERT(s2toZAlignedSO3(v3));
    REQUIRE_NO_ASSERT(s2_d_oplus_d_delta_0(v3));

    Eigen::MatrixX<S> m21(2, 1);
    REQUIRE_ASSERT(s2toZAlignedSO3(m21));
    REQUIRE_ASSERT(s2_d_oplus_d_delta_0(m21));

    Eigen::VectorX<S> v4(4);
    REQUIRE_ASSERT(s2toZAlignedSO3(v4));
    REQUIRE_ASSERT(s2_d_oplus_d_delta_0(v4));

    Eigen::VectorX<S> v10 = Eigen::VectorX<S>::Zero(10);
    v10(0) = 1;
    REQUIRE_NO_ASSERT(s2toZAlignedSO3(v10.segment(0, 3)));
    REQUIRE_ASSERT(s2toZAlignedSO3(v10.segment(0, 4)));
    REQUIRE_NO_ASSERT(s2_d_oplus_d_delta_0(v10.segment(0, 3)));
    REQUIRE_ASSERT(s2_d_oplus_d_delta_0(v10.segment(0, 4)));
  }
  if (csLie::Asserts::isAssertParanoidEnabled()) {
    REQUIRE_NO_ASSERT(s2toZAlignedSO3(Eigen::Vector3<S>::UnitX()));
    REQUIRE_ASSERT(s2toZAlignedSO3(Eigen::Vector3<S>::UnitX() * 2));
    REQUIRE_ASSERT(s2toZAlignedSO3(Eigen::Vector3<S>::Zero()));

    REQUIRE_NO_ASSERT(s2_d_oplus_d_delta_0(Eigen::Vector3<S>::UnitX()));
    REQUIRE_ASSERT(s2_d_oplus_d_delta_0(Eigen::Vector3<S>::UnitX() * 2));
    REQUIRE_ASSERT(s2_d_oplus_d_delta_0(Eigen::Vector3<S>::Zero()));
  }

  // s2oplus(Eigen::Matrix<S,10,1>(), Eigen::Matrix<S,2,1>()) // does not
  // compile s2oplus(Eigen::Matrix<S,3,1>(), Eigen::Matrix<S,10,1>()) // does
  // not compile
  if (csLie::Asserts::isAssertDebugEnabled()) {
    Eigen::MatrixX<S> m31(3, 1);
    m31.setZero();
    m31(0, 0) = 1;
    Eigen::MatrixX<S> m21(2, 1);
    m21.setZero();

    REQUIRE_NO_ASSERT(s2oplus(m31, m21));

    Eigen::MatrixX<S> m51(5, 1);
    m51.setZero();
    m51(0, 0) = 1;

    REQUIRE_ASSERT(s2oplus(m51, m21));
    REQUIRE_ASSERT(s2oplus(m31, m51));
  }
  if (csLie::Asserts::isAssertParanoidEnabled()) {
    Eigen::MatrixX<S> m31(3, 1);
    m31.setZero();
    m31(0, 0) = 1;
    Eigen::MatrixX<S> m21(2, 1);
    m21.setZero();
    REQUIRE_NO_ASSERT(s2oplus(m31, m21));
    REQUIRE_ASSERT(s2oplus(m31 * 2.0, m21));
    REQUIRE_ASSERT(s2oplus(m31 * 0., m21));
  }
}

} // namespace csLie
