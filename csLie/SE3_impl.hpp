#pragma once

#include <csCore/eigenUtils/IsometryMatrix_impl.hpp>
#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SE3.hpp>
#include <csLie/SO3_impl.hpp>
#include <csLie/csLieAssert.hpp>
namespace csLie {

template <class Derived> SE3<typename Derived::Scalar> se3Exp(Eigen::MatrixBase<Derived> const &delta) {
  using S = typename Derived::Scalar;
  using std::cos;
  using std::sin;
  using std::sqrt;
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SE3Tangent<S>()), "delta required to be 6x1");

  SE3<S> transf;

  Eigen::Matrix<S, 3, 1> const w2 = delta.template bottomLeftCorner<3, 1>().cwiseAbs2();
  S const theta_sq = w2.sum(); // sum of the squares
  S const theta = sqrt(theta_sq);

  if (theta < Constants<S>::epsilon) {
    S tmp5 = static_cast<S>(0.5) - theta_sq * (static_cast<S>(1.0 / 24.0) - static_cast<S>(1.0 / 720.0) * theta_sq);
    S tmp6 = delta(3, 0) * delta(4, 0) * tmp5;
    S tmp7 = static_cast<S>(1.0) - theta_sq * (static_cast<S>(1.0 / 120.0) * theta_sq + static_cast<S>(1.0 / 6.0));
    S tmp8 = delta(3, 0) * delta(5, 0) * tmp5;
    S tmp9 = delta(4, 0) * delta(5, 0) * tmp5;
    S tmp10 = delta(2, 0) * delta(4, 0) - delta(1, 0) * delta(5, 0);
    S tmp11 = delta(1, 0) * delta(3, 0) - delta(0, 0) * delta(4, 0);
    S tmp12 = delta(0, 0) * delta(5, 0) - delta(2, 0) * delta(3, 0);
    S tmp13 =
        static_cast<S>(1.0 / 6.0) - theta_sq * (static_cast<S>(1.0 / 120.0) - static_cast<S>(1.0 / 5040.0) * theta_sq);

    transf(0, 0) = 1 - (w2(2) + w2(1)) * tmp5;
    transf(1, 0) = delta(5, 0) * tmp7 + tmp6;
    transf(2, 0) = tmp8 - delta(4, 0) * tmp7;
    transf(0, 1) = tmp6 - delta(5, 0) * tmp7;
    transf(1, 1) = 1 - (w2(2) + w2(0)) * tmp5;
    transf(2, 1) = delta(3, 0) * tmp7 + tmp9;
    transf(0, 2) = delta(4, 0) * tmp7 + tmp8;
    transf(1, 2) = tmp9 - delta(3, 0) * tmp7;
    transf(2, 2) = 1 - (w2(1) + w2(0)) * tmp5;

    transf(0, 3) = (delta(4, 0) * tmp11 - delta(5, 0) * tmp12) * tmp13 + tmp10 * tmp5 + delta(0, 0);
    transf(1, 3) = (delta(5, 0) * tmp10 - delta(3, 0) * tmp11) * tmp13 + tmp12 * tmp5 + delta(1, 0);
    transf(2, 3) = (delta(3, 0) * tmp12 - delta(4, 0) * tmp10) * tmp13 + tmp11 * tmp5 + delta(2, 0);

  } else {
    S tmp5 = 1 / theta_sq;
    S tmp7 = 1 - cos(theta);
    S tmp8 = delta(3, 0) * delta(4, 0) * tmp5 * tmp7;
    S tmp9 = 1 / theta;
    S tmp10 = sin(theta);
    S tmp11 = delta(3, 0) * delta(5, 0) * tmp5 * tmp7;
    S tmp12 = delta(4, 0) * delta(5, 0) * tmp5 * tmp7;
    S tmp13 = delta(2, 0) * delta(4, 0) - delta(1, 0) * delta(5, 0);
    S tmp14 = delta(1, 0) * delta(3, 0) - delta(0, 0) * delta(4, 0);
    S tmp15 = delta(0, 0) * delta(5, 0) - delta(2, 0) * delta(3, 0);
    S tmp16 = 1 - tmp9 * tmp10;

    // variable: res
    transf(0, 0) = 1 - (w2(2) + w2(1)) * tmp5 * tmp7;
    transf(1, 0) = delta(5, 0) * tmp9 * tmp10 + tmp8;
    transf(2, 0) = tmp11 - delta(4, 0) * tmp9 * tmp10;
    transf(0, 1) = tmp8 - delta(5, 0) * tmp9 * tmp10;
    transf(1, 1) = 1 - (w2(2) + w2(0)) * tmp5 * tmp7;
    transf(2, 1) = delta(3, 0) * tmp9 * tmp10 + tmp12;
    transf(0, 2) = delta(4, 0) * tmp9 * tmp10 + tmp11;
    transf(1, 2) = tmp12 - delta(3, 0) * tmp9 * tmp10;
    transf(2, 2) = 1 - (w2(1) + w2(0)) * tmp5 * tmp7;

    transf(0, 3) = (delta(4, 0) * tmp14 - delta(5, 0) * tmp15) * tmp5 * tmp16 + tmp13 * tmp5 * tmp7 + delta(0, 0);
    transf(1, 3) = (delta(5, 0) * tmp13 - delta(3, 0) * tmp14) * tmp5 * tmp16 + tmp15 * tmp5 * tmp7 + delta(1, 0);
    transf(2, 3) = (delta(3, 0) * tmp15 - delta(4, 0) * tmp13) * tmp5 * tmp16 + tmp14 * tmp5 * tmp7 + delta(2, 0);
  }

  transf.matrix().row(3).template head<3>().setZero();
  transf(3, 3) = 1;

  return transf;
}

template <class Derived> SE3Tangent<typename Derived::Scalar> se3Log(Eigen::MatrixBase<Derived> const &transf) {
  using S = typename Derived::Scalar;
  using std::cos;
  using std::sin;
  using std::sqrt;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(transf, SE3Matrix<S>()),
                     "transf required to be 4x4");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidIsometryMatrix3D(transf), "transf is not a valid isometry matrix");

  // SO3 log
  SE3Tangent<S> ret;
  ret.template tail<3>() = so3Log(transf.template topLeftCorner<3, 3>());

  const S theta_sq = ret.template tail<3>().cwiseAbs2().sum();
  const S theta = sqrt(theta_sq);

  S a, b, c;
  if (theta < Constants<S>::epsilon) {
    a = 1 - theta_sq * (static_cast<S>(1.0 / 6) + theta_sq * static_cast<S>(1.0 / 120));
    b = static_cast<S>(0.5) - theta_sq * (static_cast<S>(1.0 / 24) - theta_sq * static_cast<S>(1.0 / 720));
    c = static_cast<S>(1.0 / 6) - theta_sq * (static_cast<S>(1.0 / 120) - theta_sq * static_cast<S>(1.0 / 5040));
  } else {
    const S inv_theta = 1 / theta;
    const S inv_theta_sq = 1 / theta_sq;
    a = inv_theta * sin(theta);
    b = inv_theta_sq * (1 - cos(theta));
    c = inv_theta_sq * (1 - a);
  }

  // TODO(sc) a bunch of magic numbers here...

  static constexpr S epsilon25{Constants<S>::epsilon * S(25)};
  S d;
  if (theta_sq < epsilon25) {
    d = static_cast<S>(1.0 / 12) + theta_sq * (static_cast<S>(1.0 / 720) + theta_sq * static_cast<S>(1.0 / 30240));
  } else if (theta_sq > 9) {
    d = (b - static_cast<S>(0.5) * a) / (b * theta_sq);
  } else {
    d = (b * static_cast<S>(0.5) - c) / a;
  }

  Eigen::Vector3<S> const wxt = ret.template tail<3>().cross(transf.template topRightCorner<3, 1>());
  Eigen::Vector3<S> const wxwxt = ret.template tail<3>().cross(wxt);

  ret.template head<3>() = transf.template topRightCorner<3, 1>() - S(0.5) * wxt + d * wxwxt;

  return ret;
}

template <class S> SE3Tangent<S> se3Log(SE3<S> const &transf) { return se3Log(transf.matrix()); }

template <class Derived> SE3Adjoint<typename Derived::Scalar> se3Adjoint(Eigen::MatrixBase<Derived> const &transf) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(transf, SE3Matrix<S>()),
                     "transf required to be 4x4");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidIsometryMatrix3D(transf), "transf is not a valid isometry matrix");

  SE3Adjoint<S> ad;
  ad.template topLeftCorner<3, 3>() = transf.template topLeftCorner<3, 3>();
  ad.template topRightCorner<3, 3>() =
      hat(transf.template topRightCorner<3, 1>()) * transf.template topLeftCorner<3, 3>();
  ad.template bottomRightCorner<3, 3>() = transf.template topLeftCorner<3, 3>();
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

  Eigen::Matrix<S, 3, 6> J{};
  J.template leftCols<3>() = rot;
  J.template rightCols<3>() = -rot * hat(point);
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
