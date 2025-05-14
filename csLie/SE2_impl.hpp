#pragma once

#include <csCore/eigenUtils/IsometryMatrix.hpp>
#include <csCore/eigenUtils/IsometryMatrix_impl.hpp>
#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csLie/SE2.hpp>
#include <csLie/SO2_impl.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

template <class Derived>
SE2<typename Derived::Scalar> se2Exp(Eigen::MatrixBase<Derived> const& delta) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SE2Alg<S>()), "delta required to be 3x1");
  SE2<S> Rt{};

  // note: to be able to use both eigen matrix and vector, it is required to use
  // always two indices
  S const theta{delta(2, 0)};
  Rt.linear() = so2Exp<S>(theta);

  S sinTh_Th{};
  S m1cosTh_Th{};
  if (std::abs(theta) < Eigen::NumTraits<S>::epsilon()) {
    sinTh_Th = 1;
    m1cosTh_Th = 0;
  } else {
    S const sinTh{Rt(1, 0)};
    S const cosTh{Rt(0, 0)};
    sinTh_Th = sinTh / theta;
    m1cosTh_Th = (1 - cosTh) / theta;
  }

  S const x{delta(0, 0)};
  S const y{delta(1, 0)};
  Rt.translation().x() = x * sinTh_Th - y * m1cosTh_Th;
  Rt.translation().y() = x * m1cosTh_Th + y * sinTh_Th;

  // bottom row
  Rt(2, 0) = 0;
  Rt(2, 1) = 0;
  Rt(2, 2) = 1;
  return Rt;
}

template <class Derived>
SE2<typename Derived::Scalar> xyYToSe2(Eigen::MatrixBase<Derived> const& xyY) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(xyY, Eigen::Matrix<S, 3, 1>()), "xyY required to be 3x1");

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

template <class Derived>
Eigen::Vector3<typename Derived::Scalar> se2ToxyY(Eigen::MatrixBase<Derived> const& transf) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(transf, SE2Matrix<S>()), "transf required to be 3x3");

  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidIsometryMatrix2D(transf), "transf is not a valid isometry matrix");

  Eigen::Vector3<S> xyY{transf(0, 2), transf(1, 2), so2Log(transf.template topLeftCorner<2, 2>())};
  return xyY;
}

template <class S>
Eigen::Vector3<S> se2ToxyY(SE2<S> const& transf) {
  return se2ToxyY(transf.matrix());
}

} // namespace csLie
