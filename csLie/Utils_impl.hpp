#pragma once

#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/Utils.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

template <class Derived>
NormalizedWithJacobian<typename Derived::Scalar> normalizeWithJacobian(Eigen::MatrixBase<Derived> const &v) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(v, Eigen::Matrix<S, 3, 1>()),
                     "normalizeWithJacobian, v required to be 3x1");

  CSLIE_ASSERT_PARANOID(std::abs(v.squaredNorm()) >= Constants<S>::epsilon,
                        "normalizeWithJacobian, input vector is zero");

  NormalizedWithJacobian<typename Derived::Scalar> out;

  Eigen::Matrix<S, 3, 1> const in2 = v.cwiseAbs2();
  S const normSq = in2.sum();
  S const norm = std::sqrt(normSq);
  S const norm_3_2 = normSq * norm;

  out.unit = v / norm;

  out.jacobian(0, 0) = 1 / norm - in2.x() / norm_3_2;
  out.jacobian(0, 1) = -v(0, 0) * v(1, 0) / norm_3_2;
  out.jacobian(0, 2) = -v(0, 0) * v(2, 0) / norm_3_2;

  out.jacobian(1, 0) = out.jacobian(0, 1);
  out.jacobian(1, 1) = 1 / norm - in2.y() / norm_3_2;
  out.jacobian(1, 2) = -v(1, 0) * v(2, 0) / norm_3_2;

  out.jacobian(2, 0) = out.jacobian(0, 2);
  out.jacobian(2, 1) = out.jacobian(1, 2);
  out.jacobian(2, 2) = 1 / norm - in2.z() / norm_3_2;

  return out;
}
} // namespace csLie
