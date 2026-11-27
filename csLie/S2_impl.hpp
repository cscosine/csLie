#pragma once

#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/S2.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

template <class Derived> S2<typename Derived::Scalar> s2Exp(Eigen::MatrixBase<Derived> const &delta) {
  using S = typename Derived::Scalar;
  using std::cos;
  using std::sin;
  using std::sqrt;
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, S2Tangent<S>()), "delta required to be 2x1");

  // (cos(theta), sin(theta) / theta delta), theta = |delta|: sin(theta) / theta has no cancellation, no series needed.
  // delta read once (an expression, e.g. Random(), gives new values at each read)
  S const d0 = delta(0, 0);
  S const d1 = delta(1, 0);
  S const theta = sqrt(d0 * d0 + d1 * d1);
  S const sinc = theta > S(0) ? sin(theta) / theta : S(1);
  return S2<S>(cos(theta), sinc * d0, sinc * d1);
}

template <class Derived> SO3<typename Derived::Scalar> s2Frame(Eigen::MatrixBase<Derived> const &v) {
  using S = typename Derived::Scalar;
  using std::abs;
  using std::hypot;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(v, S2<S>()), "v required to be 3x1");
  CSLIE_ASSERT_PARANOID(abs(v.squaredNorm() - S(1)) < Constants<S>::dummy_precision, "v is not a unit vector");

  // [v, b1, b2] with b2 = (a x v) / |a x v| and b1 = b2 x v, a the pivot axis: x while x^2 <= 1/2, else z. Two charts,
  // |a x v| >= 1/sqrt(2) in both: no trigonometry, well conditioned everywhere; the frame jumps at x^2 = 1/2
  S const x = v(0, 0);
  S const y = v(1, 0);
  S const z = v(2, 0);
  SO3<S> R;
  R.col(0) << x, y, z;
  if (x * x <= S(0.5)) {
    // a = x: a x v = (0, -z, y)
    S const r = hypot(y, z);
    R.col(1) << -r, x * y / r, x * z / r;
    R.col(2) << S(0), -z / r, y / r;
  } else {
    // a = z: a x v = (-y, x, 0)
    S const r = hypot(x, y);
    R.col(1) << x * z / r, y * z / r, -r;
    R.col(2) << -y / r, x / r, S(0);
  }
  return R;
}

template <class Derived1, class Derived2>
S2<typename Derived1::Scalar> s2Oplus(Eigen::MatrixBase<Derived1> const &v, Eigen::MatrixBase<Derived2> const &delta) {
  using S = typename Derived1::Scalar;
  static_assert(std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value,
                "s2Oplus needs equal scalar types");
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, S2Tangent<S>()), "delta required to be 2x1");
  return S2<S>(s2Frame(v) * s2Exp(delta.eval()));
}

template <class Derived1, class Derived2>
S2Tangent<typename Derived1::Scalar> s2OplusInverse(Eigen::MatrixBase<Derived1> const &v,
                                                    Eigen::MatrixBase<Derived2> const &w) {
  using S = typename Derived1::Scalar;
  using std::atan2;
  using std::hypot;
  static_assert(std::is_same<typename Derived1::Scalar, typename Derived2::Scalar>::value,
                "s2OplusInverse needs equal scalar types");
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(w, S2<S>()), "w required to be 3x1");

  // w in the frame of v: (cos(theta), sin(theta) d) with d the unit direction of delta. The tangent part from w - v
  // (exact for w near v, and b1, b2 orthogonal to v), accurate to a few epsilon of |delta| at small angles
  S2<S> const vv = v; // read once
  S2<S> const ww = w;
  SO3<S> const R = s2Frame(vv);
  S2<S> const wv = ww - vv;
  S const c = vv.dot(ww);
  S const t1 = R.col(1).dot(wv);
  S const t2 = R.col(2).dot(wv);
  S const sinTheta = hypot(t1, t2);
  if (!(sinTheta > S(0))) {
    // w = v, or w = -v: every direction is at the angle pi, the one of b1 (the cut locus: the direction jumps near it)
    return c > S(0) ? S2Tangent<S>::Zero() : S2Tangent<S>(Constants<S>::pi, S(0));
  }
  S const k = atan2(sinTheta, c) / sinTheta;
  return S2Tangent<S>(k * t1, k * t2);
}

template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 2> s2OplusJacobianAtZero(Eigen::MatrixBase<Derived> const &v) {
  return s2Frame(v).template rightCols<2>();
}

} // namespace csLie
