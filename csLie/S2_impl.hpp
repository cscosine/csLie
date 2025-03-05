#pragma once

#include <csLie/S2.hpp>
#include <csLie/Constants.hpp>

#include <libassert/assert.hpp>

namespace csLie {

  /// @brief a constant representing epsilon^1/4, used as boundary to evaluate s2 exp with taylor approximation
  /// @tparam S the scalar type
  template <class S> static const S s2TaylorBound{ std::sqrt(std::sqrt(Constants<S>::epsilon)) };

  template <class S> SO3<S> s2toZAlignedSO3(S2<S> const& v) {
    DEBUG_ASSERT( std::abs(v.squaredNorm() - 1.0) < Constants<S>::dummy_precision, "s2toZAlignedSO3 input vector is not unitary");

    SO3<S> R{};
    // x direction aligned with v
    R.col(0) = v;
    // right cols are equivalent to the jacobian
    R.template rightCols<2>() = s2_d_OPlus_d_delta_0(v);

    return R;
  }

  template <class S> S2<S> s2Exp(const S2Alg<S>& delta) {
    S const norm2 = delta.squaredNorm();

    // TODO(sc) verify this logic and the taylor boundary
    S cos_norm{}, sinc_norm{};
    if (norm2 > s2TaylorBound<S>) {
      S const x = std::sqrt(norm2);
      cos_norm = std::cos(x);
      sinc_norm = std::sin(x) / x;
    }
    else {
      // TODO(sc) move to a static const expr 
      static const std::array<S, 7>  inv = { S(1 / 3.), S(1 / 4.), S(1 / 5.), S(1 / 6.), S(1 / 7.), S(1 / 8.), S(1 / 9.) };
      S cosi = 1.;
      S sinci = 1.;
      S term = -1 / 2. * norm2;
      for (int i = 0; i < 3; ++i) {
        cosi += term;
        term *= inv[2 * i];
        sinci += term;
        term *= -inv[2 * i + 1] * norm2;
      }
      cos_norm = cosi;
      sinc_norm = sinci;
    }

    S2<S> ret;
    ret(0) = cos_norm;
    ret(1) = delta(0) * sinc_norm;
    ret(2) = delta(1) * sinc_norm;

    return ret;

  }

  template <class S> S2<S> s2OPlus(S2<S> const& v, S2Alg<S> const& delta) {
    SO3<S> const R = s2toZAlignedSO3(v);
    S2<S> const vDelta = s2Exp(delta);
    S2<S> const ret = R * vDelta;
    return ret;
  }

  template <class S> S2OPlusJac<S> s2_d_OPlus_d_delta_0(S2<S> const& v) {
    // TODO(sc) here we assume atan2(0,0) = 0. This is protected with a test (testAtan2). In case the test fail on some platform, we need to handle the logic here
    S const alpha = std::atan2(v.z(), v.y());
    S const r = std::sqrt(v.z() * v.z() + v.y() * v.y());
    S const c = cos(alpha);
    S const s = sin(alpha);

    S2OPlusJac<S> J{};
    J.col(0) = Eigen::Matrix<S, 3, 1>(-r, v.x() * c, v.x() * s);
    J.col(1) = Eigen::Matrix<S, 3, 1>(0, -s, c);
    return J;
  }
}