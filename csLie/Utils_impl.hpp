#pragma once

#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/Utils.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

template <class Derived>
UnitVector3dWithJacobians<typename Derived::Scalar>
normalize3dVectorWithJacobians(Eigen::MatrixBase<Derived> const &v) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(v, Eigen::Matrix<S, 3, 1>()),
                     "normalize3dVectorWithJacobians, v required to be 3x1");

  CSLIE_ASSERT_PARANOID(std::abs(v.squaredNorm()) >= Constants<S>::epsilon,
                        "normalize3dVectorWithJacobians, input vector is zero");

  UnitVector3dWithJacobians<typename Derived::Scalar> out;

  Eigen::Matrix<S, 3, 1> const in2 = v.cwiseAbs2();
  S const normSq = in2.sum();
  S const norm = std::sqrt(normSq);
  S const norm_3_2 = normSq * norm;

  out.vo = v / norm;

  out.d_vo_d_v(0, 0) = 1 / norm - in2.x() / norm_3_2;
  out.d_vo_d_v(0, 1) = -v(0, 0) * v(1, 0) / norm_3_2;
  out.d_vo_d_v(0, 2) = -v(0, 0) * v(2, 0) / norm_3_2;

  out.d_vo_d_v(1, 0) = out.d_vo_d_v(0, 1);
  out.d_vo_d_v(1, 1) = 1 / norm - in2.y() / norm_3_2;
  out.d_vo_d_v(1, 2) = -v(1, 0) * v(2, 0) / norm_3_2;

  out.d_vo_d_v(2, 0) = out.d_vo_d_v(0, 2);
  out.d_vo_d_v(2, 1) = out.d_vo_d_v(1, 2);
  out.d_vo_d_v(2, 2) = 1 / norm - in2.z() / norm_3_2;

  return out;
}
} // namespace csLie
