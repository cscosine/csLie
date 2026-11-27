#pragma once

#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csCore/eigenUtils/RotationMatrix_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SO3.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

namespace _private {

/// @brief the coefficients of the SO3 exp, R = I + A W + B W^2 (W = skewSymm(omega), theta = |omega|)
template <class S> struct SO3ExpCoeffs {
  S A; ///< sin(theta) / theta
  S B; ///< (1 - cos(theta)) / theta^2
};

/// @brief the coefficients of the SO3 exp at the angle theta >= 0, from the half angle: with
/// h = sin(theta / 2) / (theta / 2) and c = cos(theta / 2), A = h c and B = h^2 / 2. Accurate to a few epsilon at every
/// angle, with no series and no switch: sin(x) / x has no cancellation for any x > 0 (B as (1 - cos(theta)) / theta^2
/// loses all its digits below theta ~ sqrt(epsilon)); at theta = 0, h = 1.
template <class S> SO3ExpCoeffs<S> so3ExpCoeffs(S const theta) {
  using std::cos;
  using std::sin;
  S const half = theta / S(2);
  S const h = half > S(0) ? sin(half) / half : S(1);
  return {h * cos(half), S(0.5) * h * h};
}

/// @brief the coefficients of the SO3 left jacobian and of its inverse (W = skewSymm(omega), theta = |omega|):
/// J = I + B W + C W^2 = A I + B W + C omega omega^T, J^-1 = I - W / 2 + D W^2 = K I - W / 2 + D omega omega^T
template <class S> struct SO3JCoeffs {
  S A; ///< sin(theta) / theta
  S B; ///< (1 - cos(theta)) / theta^2
  S C; ///< (theta - sin(theta)) / theta^3
  S D; ///< (1 - A / (2 B)) / theta^2
  S K; ///< A / (2 B) = (theta / 2) cot(theta / 2) = 1 - D theta^2
};

/// @brief the number of terms of the series of C(x) = (x - sin(x)) / x^3 = sum_k (-x^2)^k / (2k + 3)! that makes its
/// truncation error below epsilon / 8 of C at x <= 2 (C(2) > 1/8): the first omitted term 4^n / (2n + 3)! below
/// epsilon / 64. A bound, computed in S at compile time
template <class S> constexpr int so3CSeriesTerms() {
  S term = S(1) / S(6); // the first omitted term, 4^n / (2n + 3)!, at n = 0
  int n = 0;
  while (term > Constants<S>::epsilon / S(64)) {
    ++n;
    term *= S(4) / S((2 * n + 2) * (2 * n + 3));
  }
  return n;
}

/// @brief C(x) = (x - sin(x)) / x^3 at x >= 0, given sin(x): accurate to a few epsilon. The closed form loses digits
/// to the cancellation of x - sin(x) below x ~ 2: there the series, with as many terms as the precision of S needs
template <class S> S so3CoeffC(S const x, S const sinX) {
  if (x > S(2)) {
    return (x - sinX) / (x * x * x);
  }
  constexpr int terms = so3CSeriesTerms<S>();
  // Horner: C = (1 - x^2 / (4 5) (1 - x^2 / (6 7) (1 - ...))) / 6
  S const x2 = x * x;
  S r = S(1);
  for (int k = terms - 1; k >= 1; --k) {
    r = S(1) - x2 * r / S((2 * k + 2) * (2 * k + 3));
  }
  return r / S(6);
}

/// @brief the coefficients of the SO3 left jacobian and of its inverse at the angle theta >= 0, accurate to a few
/// epsilon at every angle (J^-1 is singular at theta = 2 pi). From the half angle x = theta / 2, with no cancellation:
/// A and B as so3ExpCoeffs, C(theta) = (C(x) + A(x) B(x)) / 4, D(theta) = (B(x) - C(x)) / (4 A(x)), K = cos(x) / A(x)
/// (C(x) ~ 1/6 against A(x) B(x) ~ 1/2 and B(x) ~ 1/2; beyond theta = 2 pi, B(x) - C(x) cancels mildly, by 2.5 at
/// theta = 10). C and D in closed form, (theta - sin(theta)) / theta^3 and
/// (1 - A / (2 B)) / theta^2, lose all their digits at small angles.
template <class S> SO3JCoeffs<S> so3JCoeffs(S const theta) {
  using std::cos;
  using std::sin;
  S const half = theta / S(2);
  S const quarter = theta / S(4);
  S const sinHalf = sin(half);
  S const cosHalf = cos(half);
  S const h = half > S(0) ? sinHalf / half : S(1);             // A(theta / 2)
  S const hq = quarter > S(0) ? sin(quarter) / quarter : S(1); // A(theta / 4)
  S const bHalf = S(0.5) * hq * hq;                            // B(theta / 2)
  S const cHalf = so3CoeffC(half, sinHalf);                    // C(theta / 2)
  return {h * cosHalf, S(0.5) * h * h, (cHalf + h * bHalf) / S(4), (bHalf - cHalf) / (S(4) * h), cosHalf / h};
}

} // namespace _private

template <class Derived> SO3<typename Derived::Scalar> so3Exp(Eigen::MatrixBase<Derived> const &delta) {
  using S = typename Derived::Scalar;
  using std::sqrt;
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SO3Alg<S>()), "delta required to be 3x1");

  SO3Alg<S> const w2 = delta.cwiseAbs2();
  S const theta_sq = w2.sum(); // sum of the squares
  // below sqrt(min) theta_sq underflows to 0 or a subnormal: theta is then 0 or inexact, harmless, the coefficients
  // are 1 and 1/2 up to theta^2
  _private::SO3ExpCoeffs<S> const k = _private::so3ExpCoeffs(sqrt(theta_sq));
  S const a = k.A;
  S const b = k.B;

  SO3<S> R;

  R(0, 0) = S(1) - b * (w2(1) + w2(2));
  R(1, 1) = S(1) - b * (w2(0) + w2(2));
  R(2, 2) = S(1) - b * (w2(0) + w2(1));

  const S Bab = b * delta(0, 0) * delta(1, 0);
  const S Bac = b * delta(0, 0) * delta(2, 0);
  const S Bbc = b * delta(1, 0) * delta(2, 0);

  R(0, 1) = Bab - a * delta(2, 0);
  R(1, 0) = Bab + a * delta(2, 0);
  R(0, 2) = Bac + a * delta(1, 0);
  R(2, 0) = Bac - a * delta(1, 0);
  R(1, 2) = Bbc - a * delta(0, 0);
  R(2, 1) = Bbc + a * delta(0, 0);

  return R;
}

template <class Derived> SO3J<typename Derived::Scalar> so3J(Eigen::MatrixBase<Derived> const &delta) {
  using S = typename Derived::Scalar;
  using std::sqrt;
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SO3Alg<S>()), "delta required to be 3x1");

  // J = A I + B W + C omega omega^T
  S const x = delta(0, 0);
  S const y = delta(1, 0);
  S const z = delta(2, 0);
  _private::SO3JCoeffs<S> const k = _private::so3JCoeffs(sqrt(x * x + y * y + z * z));

  S const Cxy = k.C * x * y;
  S const Cxz = k.C * x * z;
  S const Cyz = k.C * y * z;

  SO3J<S> ret;
  ret(0, 0) = k.A + k.C * x * x;
  ret(1, 1) = k.A + k.C * y * y;
  ret(2, 2) = k.A + k.C * z * z;
  ret(0, 1) = Cxy - k.B * z;
  ret(1, 0) = Cxy + k.B * z;
  ret(0, 2) = Cxz + k.B * y;
  ret(2, 0) = Cxz - k.B * y;
  ret(1, 2) = Cyz - k.B * x;
  ret(2, 1) = Cyz + k.B * x;
  return ret;
}

template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 1> so3Log(Eigen::MatrixBase<Derived> const &rot) {
  using S = typename Derived::Scalar;
  using std::atan2;
  using std::sqrt;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(rot, SO3<S>()), "rot required to be 3x3");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidRotationMatrix3D(rot), "rot is not a valid rotation matrix");

  // Rotation matrix to a quaternion (w, v), up to a positive scale, with Shepperd's method: the square root is taken
  // of the largest of 1 + trace and the diagonal terms, so nothing is divided by a value near zero. Then
  // theta = 2 atan2(|v|, w) and omega = theta v / |v|, accurate at every angle (no acos, no series, no switch).
  // As GTSAM's SO3::Logmap (Hertzberg et al., Information Fusion, 2011).
  S const tr = rot(0, 0) + rot(1, 1) + rot(2, 2);
  S w;
  SO3Alg<S> v;
  if (tr > S(0)) {
    // theta < 2 pi / 3: (1 + tr, vee(R - R^T)) = 4 cos(theta / 2) (w, v)
    w = S(1) + tr;
    v << rot(2, 1) - rot(1, 2), rot(0, 2) - rot(2, 0), rot(1, 0) - rot(0, 1);
  } else {
    // pivot on the largest diagonal term i: 1 + 2 R(i, i) - tr = 4 v_i^2 (up to the scale)
    int i = 0;
    if (rot(1, 1) > rot(i, i)) {
      i = 1;
    }
    if (rot(2, 2) > rot(i, i)) {
      i = 2;
    }
    int const j = (i + 1) % 3;
    int const k = (i + 2) % 3;
    S const r = sqrt(S(1) + S(2) * rot(i, i) - tr);
    S const invR = S(1) / r;
    v(i) = r;
    v(j) = (rot(j, i) + rot(i, j)) * invR;
    v(k) = (rot(k, i) + rot(i, k)) * invR;
    w = (rot(k, j) - rot(j, k)) * invR;
    // theta in [0, pi]; at exactly pi (w = 0) the sign is the pivot's: deterministic
    if (w < S(0)) {
      w = -w;
      v = -v;
    }
  }

  S const n = v.norm();
  if (!(n > S(0))) {
    // the identity (the other branch has n >= r > 0)
    return SO3Alg<S>::Zero();
  }
  S const theta = S(2) * atan2(n, w);
  return (theta / n) * v;
}

template <class Derived> SO3Jinv<typename Derived::Scalar> so3Jinv(Eigen::MatrixBase<Derived> const &delta) {
  using S = typename Derived::Scalar;
  using std::sqrt;
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SO3Alg<S>()), "delta required to be 3x1");

  // J^-1 = K I - W / 2 + D omega omega^T
  S const x = delta(0, 0);
  S const y = delta(1, 0);
  S const z = delta(2, 0);
  _private::SO3JCoeffs<S> const k = _private::so3JCoeffs(sqrt(x * x + y * y + z * z));

  S const Dxy = k.D * x * y;
  S const Dxz = k.D * x * z;
  S const Dyz = k.D * y * z;
  S const hx = x / S(2);
  S const hy = y / S(2);
  S const hz = z / S(2);

  SO3Jinv<S> ret;
  ret(0, 0) = k.K + k.D * x * x;
  ret(1, 1) = k.K + k.D * y * y;
  ret(2, 2) = k.K + k.D * z * z;
  ret(0, 1) = Dxy + hz;
  ret(1, 0) = Dxy - hz;
  ret(0, 2) = Dxz - hy;
  ret(2, 0) = Dxz + hy;
  ret(1, 2) = Dyz + hx;
  ret(2, 1) = Dyz - hx;
  return ret;
}

// reference for RPY conversion
// https://eecs.qmul.ac.uk/~gslabaugh/publications/euler.pdf
template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 1> so3ToRPY(Eigen::MatrixBase<Derived> const &rot) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(rot, SO3<S>()), "rot required to be 3x3");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidRotationMatrix3D(rot), "rot is not a valid rotation matrix");

  const S sin_pitch{rot(2, 0)};
  if (std::abs(std::abs(sin_pitch) - static_cast<S>(1.0)) > Constants<S>::epsilon) {
    const S pitch{-std::asin(rot(2, 0))};
    const S cos_pitch{std::cos(pitch)}; // can we extract it from the matrix?

    const S roll{std::atan2(rot(2, 1) / cos_pitch, rot(2, 2) / cos_pitch)};
    const S yaw{std::atan2(rot(1, 0) / cos_pitch, rot(0, 0) / cos_pitch)};
    return Eigen::Matrix<S, 3, 1>(roll, pitch, yaw);
  } else {
    if (sin_pitch < 0) { // ~= -1
      const S pitch{Constants<S>::pi / static_cast<S>(2.0)};
      const S yaw = 0;
      const S roll = std::atan2(rot(0, 1), rot(0, 2));
      return Eigen::Matrix<S, 3, 1>(roll, pitch, yaw);
    } else { // ~= 1
      const S pitch{-Constants<S>::pi / static_cast<S>(2.0)};
      const S yaw = 0;
      const S roll = std::atan2(-rot(0, 1), -rot(0, 2));
      return Eigen::Matrix<S, 3, 1>(roll, pitch, yaw);
    }
  }
}

template <class Derived> SO3<typename Derived::Scalar> rpyToSo3(Eigen::MatrixBase<Derived> const &rpy) {
  using S = typename Derived::Scalar;
  using std::cos;
  using std::sin;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(rpy, Eigen::Matrix<S, 3, 1>()),
                     "rpy required to be 3x1");

  SO3<S> R;

  S const roll_x = rpy(0, 0);
  S const pitch_y = rpy(1, 0);
  S const yaw_z = rpy(2, 0);

  S const c_a = cos(yaw_z);
  S const s_a = sin(yaw_z);
  S const c_b = cos(pitch_y);
  S const s_b = sin(pitch_y);
  S const c_c = cos(roll_x);
  S const s_c = sin(roll_x);

  R(0, 0) = c_a * c_b;
  R(1, 0) = s_a * c_b;
  R(2, 0) = -s_b;

  R(0, 1) = c_a * s_b * s_c - s_a * c_c;
  R(1, 1) = s_a * s_b * s_c + c_a * c_c;
  R(2, 1) = c_b * s_c;

  R(0, 2) = c_a * s_b * c_c + s_a * s_c;
  R(1, 2) = s_a * s_b * c_c - c_a * s_c;
  R(2, 2) = c_b * c_c;

  return R;
}

template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymm(Eigen::MatrixBase<Derived> const &v) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(v, Eigen::Matrix<S, 3, 1>()),
                     "v required to be 3x1");

  Eigen::Matrix<S, 3, 3> R;

  R << 0, -v(2, 0), v(1, 0), v(2, 0), 0, -v(0, 0), -v(1, 0), v(0, 0), 0;

  return R;
}

} // namespace csLie
