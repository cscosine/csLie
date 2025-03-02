#include <catch2/catch_test_macros.hpp>

#include <csLie/Types.hpp>
#include <csLie/S2_impl.hpp>

#include <iostream>
#include <Eigen/Dense>

static const int N_RANDOM = 1000;

// TODO(sc) template test on float / double

TEST_CASE("S2 test - assumption", "[S2]") {
  // in s2toZAlignedSO3 we assume that atan2(0,0) = 0. This might be platform dependent. 
  // In case this test fail, adapt the logic in s2toZAlignedSO3 to handle the edge case
  REQUIRE(std::atan2(0, 0) == 0.0);
}

TEST_CASE("S2 test - base", "[S2]") {

  for (int i = 0; i < N_RANDOM; i++) {
    csLie::S2d v;
    v.setRandom().normalize();

    csLie::SO3d Rv = csLie::s2toZAlignedSO3(v);
    // TODO(sc) check is a proper rot matrix (add a function common for that)
    // TODO(sc) check Rx = v
    REQUIRE(std::abs(Rv.determinant() - 1) < Eigen::NumTraits<double>::dummy_precision());

    csLie::S2Algd delta;
    delta.setRandom();
    csLie::S2d vDelta= csLie::s2Exp(delta);
    // TODO check edge cases (taylor expansion for small numbers)
    REQUIRE(std::abs(vDelta.norm() - 1) < Eigen::NumTraits<double>::dummy_precision());

    // oplus(v,delta)
    csLie::S2d v_oplus_expected{Rv * vDelta};
    REQUIRE(std::abs(v_oplus_expected.norm() - 1) < Eigen::NumTraits<double>::dummy_precision());

    // oplus computed with function, has to be equivalent
    csLie::S2d v_oplus{csLie::s2OPlus(v, delta) };
    REQUIRE((v_oplus - v_oplus_expected).norm() < Eigen::NumTraits<double>::dummy_precision());

    csLie::S2OPlusJacd J = csLie::s2_d_OPlus_d_delta_0(v);
    REQUIRE((J - Rv.rightCols<2>()).squaredNorm() == 0);

  }

}

/*
TEST_CASE("Normal Estimation - point plane", "[NormPPl]") {
  Eigen::Isometry3d Rt = Eigen::Isometry3d::Identity();
  Rt.translation().setRandom();
  Rt.linear().col(0).setRandom().normalize();
  Rt.linear().col(1).setRandom().normalize();
  Rt.linear().col(2) = Rt.linear().col(0).cross(Rt.linear().col(1)).normalized();
  Rt.linear().col(1) = -Rt.linear().col(0).cross(Rt.linear().col(2)).normalized();

  std::cout << "Rt " << std::endl << Rt.matrix() << std::endl;

  int N = 4;
  Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, N);

  Eigen::Vector3d pointCenter = Rt.translation();
  Eigen::Vector3d normalGT = Rt.linear().col(2);

  for (int i = 0; i < N; i++) {
    Eigen::Vector3d pr = Eigen::Vector3d::Random();
    pr.z() = 0;
    points.col(i) = Rt * pr;
    double ep = (pointCenter - points.col(i)).dot(normalGT);
    assert(std::abs(ep) < Eigen::NumTraits<double>::dummy_precision());
  }

  Eigen::Matrix<double, 3, Eigen::Dynamic> pointsNoised = points;
  Eigen::Vector3d estNormal = normalGT;
  bool addNoise = true;
  if (addNoise) {
    pointsNoised += Eigen::Matrix<double, 3, Eigen::Dynamic>::Random(3, N) * 0.02;
    estNormal += Eigen::Vector3d::Random() * 0.01;
  }

  double chi2 = 0;
  for (int it = 0; it < 10; it++) {
    Eigen::Matrix2d H = Eigen::Matrix2d::Zero();
    Eigen::Vector2d b = Eigen::Vector2d::Zero();
    chi2 = 0;
    for (int i = 0; i < N; i++) {
      Eigen::Vector3d pd = (pointCenter - pointsNoised.col(i));
      double e = pd.dot(estNormal);
      Eigen::Matrix<double, 1, 2> d_e_d_delta = pd.transpose() * csLie::d_Rv_expd_d_expd0(csLie::S2d{estNormal});

      double c2 = e * e;
      double w = 1;

      H += d_e_d_delta.transpose() * w * d_e_d_delta;
      b += d_e_d_delta.transpose() * w * e;
      chi2 += c2;
    }
    std::cout << "it " << it << ", chi2 " << chi2 << std::endl;
    Eigen::LDLT<Eigen::Matrix2d> ldlt;
    ldlt.compute(H);
    Eigen::Vector2d x = ldlt.solve(-b);
    estNormal = csLie::s2toRotationMatrix(csLie::S2d{estNormal})* csLie::exp(csLie::S2Algd{x});
    estNormal.normalize();
  }

  REQUIRE(chi2 < 1e-3);

  double err = estNormal.dot(normalGT);
  std::cout << "err wrt GT (should be close to 1)" << err << std::endl;
  REQUIRE(std::abs(1.0 - err) < 1e-3);

}
*/