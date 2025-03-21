#include <catch2/catch_all.hpp>

#include <csLie/Types.hpp>
#include <csLie/S2_impl.hpp>
#include <csLie/Constants.hpp>

#include <csLie/tests/TestUtils.hpp>

#include <Eigen/Dense> // for determinant

namespace csLie {
  

  TEMPLATE_TEST_CASE("S2 test - verify assumption atan2(0,0) = 0", "[S2]", float32_t, float64_t) {
    // in s2toZAlignedSO3 we assume that atan2(0,0) = 0. This might be platform dependent. 
    // In case this test fail, adapt the logic in s2toZAlignedSO3 to handle this edge case
    REQUIRE(std::atan2<TestType>(0, 0) == 0.0);
  }

  TEMPLATE_TEST_CASE("S2 test - base", "[S2]", float32_t, float64_t) {
    using S = TestType;
    S const PRECISION = Constants<S>::dummy_precision;

    // create a unit vector v
    S2<S> const v{ S2<S>::Random().normalized() };

    // verify that v is a unit vector (test self-sanity)
    REQUIRE_THAT(v.norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));

    // create rotation matrix from v
    SO3<S> const Rv = s2toZAlignedSO3<S>(v);
    // check is a valid rot matrix
    checkValidRotationMatrix<S>(Rv);
    // and its first column is v
    REQUIRE((Rv.col(0) - v).isZero());

    // create a random s2 algebra vector delta
    S2Alg<S> const delta{ S2Alg<S>::Random() };
    
    // create a unit vector vd as exp of delta
    S2<S> const vd = s2Exp<S>(delta);
    // TODO(sc) check edge cases (taylor expansion for small numbers)
    REQUIRE_THAT(vd.norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));

    // compute oplus(v,delta) explicitely in voe
    S2<S> const voe{ Rv * vd };
    // verify voe is a unit vector
    REQUIRE_THAT(voe.norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));

    // compute oplus(v,delta) with s2OPlus in vo
    S2<S> const vo{ s2OPlus<S>(v, delta) };
    // verify vo is a unit vector
    REQUIRE_THAT(vo.norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));
    // verify vo is equivalent to voe
    REQUIRE((vo - voe).isZero());

    // compute the jacobian of oplus(v,delta) wrt delta, evalauted in delta = 0
    S2OPlusJac<S> const J = s2_d_OPlus_d_delta_0<S>(v);
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
    Eigen::Transform<S, 3, Eigen::Isometry> Rt{ Eigen::Transform<S, 3, Eigen::Isometry> ::Identity()};
    Rt.translation().setRandom();
    Rt.linear() = Eigen::Matrix<S,3,3>(
        Eigen::AngleAxis<S>(Eigen::Matrix<S,1,1>::Random()(0), Eigen::Vector3<S>::UnitX())
      * Eigen::AngleAxis<S>(Eigen::Matrix<S, 1, 1>::Random()(0), Eigen::Vector3<S>::UnitY())
      * Eigen::AngleAxis<S>(Eigen::Matrix<S, 1, 1>::Random()(0), Eigen::Vector3<S>::UnitZ()))      ;

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
      // i.e., a vector from the center of the reference system to the point is orthogonal to the plane normal
      S const ep = (pointCenter - points.col(i)).dot(normalGT);
      REQUIRE_THAT(ep, Catch::Matchers::WithinAbs(0.0, PRECISION));
    }

    Eigen::Matrix<S, 3, Eigen::Dynamic> pointsNoised = points;
    if(pointsNoisedFlag) {
      // eventually, add noise to points
      pointsNoised += Eigen::Matrix<S, 3, Eigen::Dynamic>::Random(3, N) * 0.02;
    }

    Eigen::Vector3<S> estNormal = normalGT;
    if(perturbedInitialization) {
      // eventually, add noise to initial estimation
      estNormal = s2OPlus<S>(estNormal, Eigen::Vector2<S>::Random() * 0.01);
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
    // the point plane error e is  e = p dot pi = [p 1] dot [n d] = p^t n + d = p^t n -c^t n = (p - c)^t * n
    // the derivative d_e_d_delta evaluated in zero, computed on (p - c)^t * oplus(n, delta), is (p - c)^t * s2_d_OPlus_d_delta_0(n)
    
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
        Eigen::Matrix<S, 1, 2> d_e_d_delta = pd.transpose() * s2_d_OPlus_d_delta_0(estNormal);

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
      estNormal = s2OPlus<S>(estNormal, x );
      estNormal.normalize();
    }

    // final checks
    if (!pointsNoisedFlag) {
      // if no noise, expect to converge exactly at zero
      REQUIRE_THAT(chi2, Catch::Matchers::WithinAbs(0.0, PRECISION));
      REQUIRE((estNormal - normalGT).isZero());
      // if no noise, estimate of normal is exact, even when the start is perturbed
      S const dotFinal = estNormal.dot(normalGT);
      REQUIRE_THAT(dotFinal, Catch::Matchers::WithinAbs(1.0, PRECISION));
    }
    else {
      // if noise we can just ask that the chi2 is decreasing
      REQUIRE(chi2 < initial_chi2);
    }

    if (perturbedInitialization) {
      // if the init is perturbed, we expect that the final estimation is closer to GT
      S const dotInitial = estNormalInitial.dot(normalGT);
      S const dotFinal = estNormal.dot(normalGT);
      REQUIRE(dotFinal > dotInitial);
    }
  
  }
}