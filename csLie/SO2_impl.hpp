#pragma once

#include <csLie/SO2.hpp>

namespace csLie {

  template <class DerivedIn, class DerivedOut> void so2Exp(Eigen::MatrixBase<DerivedIn> const& delta, Eigen::MatrixBase<DerivedOut>& res)
  {
    using S = typename DerivedIn::Scalar;
    using M11 = Eigen::Matrix<S, 1, 1>;
    EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(DerivedIn, M11);
    using M22 = Eigen::Matrix<S, 2, 2>;
    EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(DerivedOut, M22);

    S const cd{ std::cos(delta(0)) };
    S const sd{ std::sin(delta(0)) };

    res << cd, -sd,
      sd, cd;

  }

  template <class S, class DerivedOut, class RetType>
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
    so2Exp(S delta, Eigen::MatrixBase<DerivedOut>& res)
  {
    so2Exp(SO2Alg<S>(delta), res);
  }

  template <class Derived> SO2<typename Derived::Scalar> so2Exp(Eigen::MatrixBase<Derived> const& delta) {
    using S = typename Derived::Scalar;
    using M11 = Eigen::Matrix<S, 1, 1>;
    EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Derived, M11);

    SO2<S> r{};
    so2Exp(delta, r);
    return r;
  }

  template <class S, class RetType>
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta) {
    return so2Exp(SO2Alg<S>(delta));
  }

}