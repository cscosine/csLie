// The left and right versions of the SO3 / SE3 jacobians and perturbations: each against its definition (central
// finite differences) and against its pair, in float32_t and float64_t. Two ways of computing the same value are
// compared within rounding, not bit for bit: the compiler may fuse a multiply and an add in one and not in the other
// (GCC on arm64 does by default).
#include <cmath>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SE3_impl.hpp>
#include <csLie/SO3_impl.hpp>
#include <csLie/Types.hpp>
#include <vector>

namespace csLie {

namespace {

/// the step and the tolerance of a central difference: step eps^(1/3), error about eps^(2/3)
template <class S> S fdStep() {
  using std::cbrt;
  return cbrt(Constants<S>::epsilon);
}
template <class S> S fdTolerance() { return S(100) * fdStep<S>() * fdStep<S>(); }

template <class S> std::vector<SO3Tangent<S>> so3Deltas() {
  return {SO3Tangent<S>(S(0), S(0), S(0)), SO3Tangent<S>(S(1e-3), S(-2e-3), S(5e-4)),
          SO3Tangent<S>(S(0.3), S(-0.7), S(0.2)), SO3Tangent<S>(S(-1.9), S(1.2), S(0.4))};
}

} // namespace

TEMPLATE_TEST_CASE("SO3 test - hat and vee", "[SO3]", float32_t, float64_t) {
  using S = TestType;
  Eigen::Vector3<S> const v(S(0.3), S(-1.2), S(2.5));
  Eigen::Vector3<S> const w(S(-0.7), S(0.4), S(1.1));
  REQUIRE((vee(hat(v)) - v).isZero(0));
  REQUIRE((hat(v) * w - v.cross(w)).isZero());
  REQUIRE((hat(v) + hat(v).transpose()).isZero(0));
}

TEMPLATE_TEST_CASE("SO3 test - left and right jacobians", "[SO3]", float32_t, float64_t) {
  using S = TestType;
  S const h = fdStep<S>();
  for (SO3Tangent<S> const &delta : so3Deltas<S>()) {
    Eigen::Matrix3<S> const Jl = so3LeftJacobian(delta);
    Eigen::Matrix3<S> const Jr = so3RightJacobian(delta);
    INFO("delta " << delta.transpose());
    // the pairs
    REQUIRE((Jr - Jl.transpose()).isZero());
    REQUIRE((Jr - so3LeftJacobian(SO3Tangent<S>(-delta))).isZero());
    REQUIRE((so3RightJacobianInverse(delta) - so3LeftJacobianInverse(delta).transpose()).isZero());
    REQUIRE((Jr * so3RightJacobianInverse(delta) - Eigen::Matrix3<S>::Identity()).isZero());
    // the definitions: exp(delta + d) = exp(Jl d) exp(delta) = exp(delta) exp(Jr d)
    SO3<S> const R = so3Exp(delta);
    Eigen::Matrix3<S> JlFd;
    Eigen::Matrix3<S> JrFd;
    for (int i = 0; i < 3; ++i) {
      SO3Tangent<S> const d = h * SO3Tangent<S>::Unit(i);
      SO3<S> const Rp = so3Exp(SO3Tangent<S>(delta + d));
      SO3<S> const Rm = so3Exp(SO3Tangent<S>(delta - d));
      JlFd.col(i) = (so3Log(Rp * R.transpose()) - so3Log(Rm * R.transpose())) / (S(2) * h);
      JrFd.col(i) = (so3Log(R.transpose() * Rp) - so3Log(R.transpose() * Rm)) / (S(2) * h);
    }
    INFO("Jl\n" << Jl << "\nfinite differences\n" << JlFd);
    REQUIRE((Jl - JlFd).cwiseAbs().maxCoeff() <= fdTolerance<S>());
    REQUIRE((Jr - JrFd).cwiseAbs().maxCoeff() <= fdTolerance<S>());
  }
}

TEMPLATE_TEST_CASE("SO3 test - oplus left and right", "[SO3]", float32_t, float64_t) {
  using S = TestType;
  SO3<S> const R = so3Exp(SO3Tangent<S>(S(0.4), S(-1.1), S(0.6)));
  for (SO3Tangent<S> const &delta : so3Deltas<S>()) {
    REQUIRE((so3OplusLeft(R, delta) - so3Exp(delta) * R).isZero());
    REQUIRE((so3OplusRight(R, delta) - R * so3Exp(delta)).isZero());
    // the left perturbation is the right one moved to the world frame: exp(R delta) R = R exp(delta)
    REQUIRE((so3OplusLeft(R, SO3Tangent<S>(R * delta)) - so3OplusRight(R, delta)).isZero());
  }
}

TEMPLATE_TEST_CASE("SE3 test - point jacobians and oplus, left and right", "[SE3]", float32_t, float64_t) {
  using S = TestType;
  S const h = fdStep<S>();
  SE3Tangent<S> xi;
  xi << S(0.5), S(-1.5), S(2), S(0.4), S(-1.1), S(0.6);
  SE3<S> const T = se3Exp(xi);
  Point3d<S> const p(S(1.2), S(-0.3), S(2.2));
  Point3d<S> const Tp = T * p;

  Eigen::Matrix<S, 3, 6> const Jl = se3PointJacobianLeft(Tp);
  Eigen::Matrix<S, 3, 6> const Jr = se3PointJacobianRight(T.linear(), p);
  Eigen::Matrix<S, 3, 6> JlFd;
  Eigen::Matrix<S, 3, 6> JrFd;
  for (int i = 0; i < 6; ++i) {
    SE3Tangent<S> const d = h * SE3Tangent<S>::Unit(i);
    SE3Tangent<S> const md = -d;
    REQUIRE((se3OplusLeft(T, d).matrix() - (se3Exp(d) * T).matrix()).isZero());
    REQUIRE((se3OplusRight(T, d).matrix() - (T * se3Exp(d)).matrix()).isZero());
    JlFd.col(i) = (se3OplusLeft(T, d) * p - se3OplusLeft(T, md) * p) / (S(2) * h);
    JrFd.col(i) = (se3OplusRight(T, d) * p - se3OplusRight(T, md) * p) / (S(2) * h);
  }
  INFO("Jl\n" << Jl << "\nfinite differences\n" << JlFd << "\nJr\n" << Jr << "\nfinite differences\n" << JrFd);
  REQUIRE((Jl - JlFd).cwiseAbs().maxCoeff() <= fdTolerance<S>() * S(10));
  REQUIRE((Jr - JrFd).cwiseAbs().maxCoeff() <= fdTolerance<S>() * S(10));
  // the pair: a right perturbation is the left one through the adjoint, Jr = Jl Ad(T)
  REQUIRE((Jr - Jl * se3Adjoint(T)).isZero());
}

TEMPLATE_TEST_CASE("SE3 test - se3PointJacobianRight, wrong sizes", "[SE3]", float32_t, float64_t) {
  using S = TestType;
  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(se3PointJacobianRight(Eigen::MatrixX<S>::Identity(3, 3), Eigen::MatrixX<S>::Zero(3, 1)));
    REQUIRE_ASSERT(se3PointJacobianRight(Eigen::MatrixX<S>::Identity(2, 3), Eigen::MatrixX<S>::Zero(3, 1)));
    REQUIRE_ASSERT(se3PointJacobianRight(Eigen::MatrixX<S>::Identity(3, 3), Eigen::MatrixX<S>::Zero(2, 1)));
  }
}

} // namespace csLie
