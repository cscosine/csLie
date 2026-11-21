#pragma once

#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csCore/eigenUtils/RotationMatrix_impl.hpp>
#include <csLie/SO2.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

template <class Derived> SO2<typename Derived::Scalar> so2Exp(Eigen::MatrixBase<Derived> const &delta) {
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SO2Alg<typename Derived::Scalar>()),
                     "delta required to be 1x1");
  // note: to be able to use both eigen matrix and vector, it is required to use
  // always two indices
  return so2Exp(delta(0, 0));
}

template <class S, class RetType> std::enable_if_t<std::is_floating_point<S>::value, RetType> so2Exp(S delta) {
  SO2<S> res{};
  S const cd{std::cos(delta)};
  S const sd{std::sin(delta)};

  res(0, 0) = cd;
  res(0, 1) = -sd;
  res(1, 0) = sd;
  res(1, 1) = cd;
  return res;
}

template <class Derived> typename Derived::Scalar so2Log(Eigen::MatrixBase<Derived> const &rot) {
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(rot, SO2<typename Derived::Scalar>()),
                     "rot required to be 2x2");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidRotationMatrix2D(rot), "rot is not a valid rotation matrix");

  // atan2 (sin,cos), average of values
  return std::atan2(rot(1, 0), rot(0, 0));
}

} // namespace csLie
