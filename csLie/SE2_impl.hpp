#pragma once

#include <csCore/eigenUtils/IsometryMatrix.hpp>
#include <csCore/eigenUtils/IsometryMatrix_impl.hpp>
#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csLie/Coefficients_impl.hpp>
#include <csLie/SE2.hpp>
#include <csLie/SO2_impl.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

template <class Derived> SE2<typename Derived::Scalar> se2Exp(Eigen::MatrixBase<Derived> const &delta) {
  using S = typename Derived::Scalar;
  using std::abs;
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SE2Alg<S>()), "delta required to be 3x1");

  // R = so2Exp(theta), t = V rho with V = [A, -B theta; B theta, A] (A, B even in theta)
  S const x{delta(0, 0)};
  S const y{delta(1, 0)};
  S const theta{delta(2, 0)};
  _private::ExpCoeffs<S> const k = _private::expCoeffs(abs(theta));
  S const bTheta = k.B * theta;

  SE2<S> Rt{};
  Rt.linear() = so2Exp<S>(theta);
  Rt.translation().x() = k.A * x - bTheta * y;
  Rt.translation().y() = bTheta * x + k.A * y;
  Rt(2, 0) = 0;
  Rt(2, 1) = 0;
  Rt(2, 2) = 1;
  return Rt;
}

template <class Derived> SE2Alg<typename Derived::Scalar> se2Log(Eigen::MatrixBase<Derived> const &transf) {
  using S = typename Derived::Scalar;
  using std::abs;
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(transf, SE2Matrix<S>()),
                     "transf required to be 3x3");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidIsometryMatrix2D(transf), "transf is not a valid isometry matrix");

  // theta = so2Log(R) in (-pi, pi], rho = V^-1 t with V^-1 = [K, theta / 2; -theta / 2, K], K = (theta / 2) cot(theta /
  // 2) (det V = 2 B): no division by a small value at any angle
  S const theta = so2Log(transf.template topLeftCorner<2, 2>());
  _private::ExpCoeffs<S> const k = _private::expCoeffs(abs(theta));
  S const kk = k.A / (S(2) * k.B);
  S const halfTheta = theta / S(2);
  S const tx = transf(0, 2);
  S const ty = transf(1, 2);
  return SE2Alg<S>(kk * tx + halfTheta * ty, -halfTheta * tx + kk * ty, theta);
}

template <class S> SE2Alg<S> se2Log(SE2<S> const &transf) { return se2Log(transf.matrix()); }

template <class Derived> SE2<typename Derived::Scalar> xyYToSe2(Eigen::MatrixBase<Derived> const &xyY) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(xyY, Eigen::Matrix<S, 3, 1>()),
                     "xyY required to be 3x1");

  SE2<S> Rt{};

  // rotation
  S const theta{xyY(2, 0)};
  Rt.linear() = so2Exp<S>(theta);

  // translation
  Rt(0, 2) = xyY(0, 0);
  Rt(1, 2) = xyY(1, 0);
  // bottom row
  Rt(2, 0) = 0;
  Rt(2, 1) = 0;
  Rt(2, 2) = 1;
  return Rt;
}

template <class Derived> Eigen::Vector3<typename Derived::Scalar> se2ToxyY(Eigen::MatrixBase<Derived> const &transf) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(transf, SE2Matrix<S>()),
                     "transf required to be 3x3");

  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidIsometryMatrix2D(transf), "transf is not a valid isometry matrix");

  Eigen::Vector3<S> xyY{transf(0, 2), transf(1, 2), so2Log(transf.template topLeftCorner<2, 2>())};
  return xyY;
}

template <class S> Eigen::Vector3<S> se2ToxyY(SE2<S> const &transf) { return se2ToxyY(transf.matrix()); }

} // namespace csLie
