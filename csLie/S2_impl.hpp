#pragma once

#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/S2.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

namespace _private {
/// @brief a constant representing epsilon^1/4, used as boundary to evaluate s2
/// exp with taylor approximation
/// @tparam S the scalar type
template <class S> static const S s2TaylorBound{std::sqrt(std::sqrt(Constants<S>::epsilon))};

} // namespace _private

template <class Derived> S2<typename Derived::Scalar> s2Exp(Eigen::MatrixBase<Derived> const &delta) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, S2Alg<S>()), "delta required to be 2x1");

  S const norm2 = delta.squaredNorm();

  // TODO(sc) verify this logic and the taylor boundary
  S cos_norm{}, sinc_norm{};
  if (norm2 > _private::s2TaylorBound<S>) {
    S const x = std::sqrt(norm2);
    cos_norm = std::cos(x);
    sinc_norm = std::sin(x) / x;
  } else {
    static constexpr size_t N{7};
    static constexpr size_t HALF_N{N / 2};
    static_assert(HALF_N == 3, "check N / 2 = 3 failed... did you changed it?");

    static constexpr std::array<S, N> inv = {S(1 / 3.), S(1 / 4.), S(1 / 5.), S(1 / 6.),
                                             S(1 / 7.), S(1 / 8.), S(1 / 9.)};
    S cosi = S(1);
    S sinci = S(1);
    S term = S(-0.5) * norm2;
    for (size_t i = 0; i < HALF_N; ++i) {
      cosi += term;
      term *= inv[2 * i];
      sinci += term;
      term *= -inv[2 * i + 1] * norm2;
    }
    cos_norm = cosi;
    sinc_norm = sinci;
  }

  S2<S> ret{};
  ret(0) = cos_norm;
  ret(1) = delta(0, 0) * sinc_norm;
  // note: to be able to use both eigen matrix and vector, it is required to use
  // always two indices accessing delta
  ret(2) = delta(1, 0) * sinc_norm;

  return ret;
}

template <class Derived> SO3<typename Derived::Scalar> s2toZAlignedSO3(Eigen::MatrixBase<Derived> const &v) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(v, S2<S>()), "v required to be 3x1");

  CSLIE_ASSERT_PARANOID(std::abs(v.squaredNorm() - S(1)) < Constants<S>::dummy_precision,
                        "s2toZAlignedSO3 input vector is not unitary");

  SO3<S> R{};
  // x direction aligned with v
  R.col(0) = v;
  // right cols are equivalent to the jacobian
  R.template rightCols<2>() = s2_d_oplus_d_delta_0(v);

  return R;
}

template <class Derived1, class Derived2>
S2<typename Derived1::Scalar> s2oplus(Eigen::MatrixBase<Derived1> const &v, Eigen::MatrixBase<Derived2> const &delta) {
  using S = typename Derived1::Scalar;
  static_assert(std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value,
                "s2oplus needs equal scalar types");
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(v, S2<S>()), "v required to be 3x1");
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, S2Alg<S>()), "delta required to be 2x1");

  CSLIE_ASSERT_PARANOID(std::abs(v.squaredNorm() - S(1)) < Constants<S>::dummy_precision,
                        "s2toZAlignedSO3 input vector is not unitary");

  SO3<S> const R = s2toZAlignedSO3(v);
  S2<S> const vDelta = s2Exp(delta);
  S2<S> const ret = R * vDelta;
  return ret;
}

template <class Derived>
s2oplusJac<typename Derived::Scalar> s2_d_oplus_d_delta_0(Eigen::MatrixBase<Derived> const &v) {
  using S = typename Derived::Scalar;
  using std::cos;
  using std::sin;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(v, S2<S>()), "v required to be 3x1");

  CSLIE_ASSERT_PARANOID(std::abs(v.squaredNorm() - S(1)) < Constants<S>::dummy_precision,
                        "s2toZAlignedSO3 input vector is not unitary");

  // note: use always two indices to allow access to v as a generic matrix
  S const x{v(0, 0)};
  S const y{v(1, 0)};
  S const z{v(2, 0)};
  // here we assume atan2(0,0) = 0. This is protected with a test (testAtan2).
  // In case the test fail on some platform, we need to handle the logic here
  S const alpha = std::atan2(z, y);
  S const r = std::sqrt(z * z + y * y);
  S const c = cos(alpha);
  S const s = sin(alpha);

  s2oplusJac<S> J{};
  J.col(0) = Eigen::Matrix<S, 3, 1>(-r, x * c, x * s);
  J.col(1) = Eigen::Matrix<S, 3, 1>(0, -s, c);
  return J;
}
} // namespace csLie
