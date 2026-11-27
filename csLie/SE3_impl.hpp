#pragma once

#include <csCore/eigenUtils/IsometryMatrix_impl.hpp>
#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csLie/Coefficients_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SE3.hpp>
#include <csLie/SO3_impl.hpp>
#include <csLie/csLieAssert.hpp>
namespace csLie {

template <class Derived> SE3<typename Derived::Scalar> se3Exp(Eigen::MatrixBase<Derived> const &delta) {
  using S = typename Derived::Scalar;
  using std::sqrt;
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SE3Tangent<S>()), "delta required to be 6x1");

  // delta = (rho, omega): R = so3Exp(omega), t = Jl(omega) rho = A rho + B omega x rho + C omega (omega . rho), the
  // coefficients of the angle |omega| accurate at every angle (Coefficients_impl.hpp)
  SE3Tangent<S> const d = delta; // read once (an expression, e.g. Random(), gives new values at each read)
  SO3Tangent<S> const rho = d.template head<3>();
  SO3Tangent<S> const omega = d.template tail<3>();
  _private::JCoeffs<S> const k = _private::jCoeffs(sqrt(omega.squaredNorm()));

  SE3<S> transf;
  transf.linear() = _private::so3ExpFromCoeffs(omega, k.A, k.B);
  transf.translation() = k.A * rho + k.B * omega.cross(rho) + (k.C * omega.dot(rho)) * omega;
  transf.matrix().row(3) << S(0), S(0), S(0), S(1);
  return transf;
}

template <class Derived> SE3Tangent<typename Derived::Scalar> se3Log(Eigen::MatrixBase<Derived> const &transf) {
  using S = typename Derived::Scalar;
  using std::sqrt;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(transf, SE3Matrix<S>()),
                     "transf required to be 4x4");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidIsometryMatrix3D(transf), "transf is not a valid isometry matrix");

  // omega = so3Log(R), rho = Jl(omega)^-1 t = K t - omega x t / 2 + D omega (omega . t); |omega| <= pi, far from the
  // singularity of Jl^-1 at 2 pi
  SE3Matrix<S> const T = transf; // read once
  SO3Tangent<S> const t = T.template topRightCorner<3, 1>();
  SO3Tangent<S> const omega = so3Log(T.template topLeftCorner<3, 3>());
  _private::JCoeffs<S> const k = _private::jCoeffs(sqrt(omega.squaredNorm()));

  SE3Tangent<S> ret;
  ret.template head<3>() = k.K * t - omega.cross(t) / S(2) + (k.D * omega.dot(t)) * omega;
  ret.template tail<3>() = omega;
  return ret;
}

template <class S> SE3Tangent<S> se3Log(SE3<S> const &transf) { return se3Log(transf.matrix()); }

template <class Derived> SE3Adjoint<typename Derived::Scalar> se3Adjoint(Eigen::MatrixBase<Derived> const &transf) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(transf, SE3Matrix<S>()),
                     "transf required to be 4x4");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidIsometryMatrix3D(transf), "transf is not a valid isometry matrix");

  SE3Matrix<S> const T = transf; // read once
  SO3<S> const R = T.template topLeftCorner<3, 3>();
  SE3Adjoint<S> ad;
  ad.template topLeftCorner<3, 3>() = R;
  ad.template topRightCorner<3, 3>() = hat(T.template topRightCorner<3, 1>()) * R;
  ad.template bottomRightCorner<3, 3>() = R;
  ad.template bottomLeftCorner<3, 3>().setZero();
  return ad;
}

template <class S> SE3Adjoint<S> se3Adjoint(SE3<S> const &transf) { return se3Adjoint(transf.matrix()); }

template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 6> se3PointJacobianLeft(Eigen::MatrixBase<Derived> const &transformedPoint) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(transformedPoint, Point3d<S>()),
                     "transformedPoint required to be 3x1");

  Eigen::Matrix<S, 3, 6> J{};
  J.template topLeftCorner<3, 3>().setIdentity();
  J.template topRightCorner<3, 3>() = -hat(transformedPoint);
  return J;
}

template <class Derived1, class Derived2>
Eigen::Matrix<typename Derived1::Scalar, 3, 6> se3PointJacobianRight(Eigen::MatrixBase<Derived1> const &rot,
                                                                     Eigen::MatrixBase<Derived2> const &point) {
  using S = typename Derived1::Scalar;
  static_assert(std::is_same<S, typename Derived2::Scalar>::value, "se3PointJacobianRight needs equal scalar types");
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(rot, SO3<S>()), "rot required to be 3x3");
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(point, Point3d<S>()), "point required to be 3x1");

  SO3<S> const R = rot; // read once
  Eigen::Matrix<S, 3, 6> J{};
  J.template leftCols<3>() = R;
  J.template rightCols<3>() = -R * hat(point);
  return J;
}

template <class S, class Derived> SE3<S> se3OplusLeft(SE3<S> const &transf, Eigen::MatrixBase<Derived> const &delta) {
  static_assert(std::is_same<S, typename Derived::Scalar>::value, "se3OplusLeft needs equal scalar types");
  return SE3<S>(se3Exp(delta) * transf);
}

template <class S, class Derived> SE3<S> se3OplusRight(SE3<S> const &transf, Eigen::MatrixBase<Derived> const &delta) {
  static_assert(std::is_same<S, typename Derived::Scalar>::value, "se3OplusRight needs equal scalar types");
  return SE3<S>(transf * se3Exp(delta));
}

} // namespace csLie
