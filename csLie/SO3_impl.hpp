#pragma once

#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csCore/eigenUtils/RotationMatrix_impl.hpp>
#include <csLie/Coefficients_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SO3.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

template <class Derived> SO3<typename Derived::Scalar> so3Exp(Eigen::MatrixBase<Derived> const &delta) {
  using S = typename Derived::Scalar;
  using std::sqrt;
  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SO3Alg<S>()), "delta required to be 3x1");

  SO3Alg<S> const w2 = delta.cwiseAbs2();
  S const theta_sq = w2.sum(); // sum of the squares
  // below sqrt(min) theta_sq underflows to 0 or a subnormal: theta is then 0 or inexact, harmless, the coefficients
  // are 1 and 1/2 up to theta^2
  _private::ExpCoeffs<S> const k = _private::expCoeffs(sqrt(theta_sq));
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
  _private::JCoeffs<S> const k = _private::jCoeffs(sqrt(x * x + y * y + z * z));

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
  _private::JCoeffs<S> const k = _private::jCoeffs(sqrt(x * x + y * y + z * z));

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
