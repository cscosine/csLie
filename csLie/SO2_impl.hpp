#pragma once

#include <csLie/SO2.hpp>
#include <csLie/SO2.hpp>
#include <csCore/eigenUtils/MatrixShape.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

  template <class Derived> 
  SO2<typename Derived::Scalar> so2Exp(Eigen::MatrixBase<Derived> const& delta) {
    using S = typename Derived::Scalar;
    using M11 = Eigen::Matrix<S, 1, 1>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, M11()), "delta required to be 1x1");
    return so2Exp(delta(0,0));
  }

  template <class S, class RetType>
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta) {
    SO2<S> res{};
    S const cd{ std::cos(delta) };
    S const sd{ std::sin(delta) };

    res(0,0) = cd;
    res(0,1) = -sd;
    res(1,0) = sd;
    res(1,1) = cd;
    return res;
  }

}