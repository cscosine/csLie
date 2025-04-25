#pragma once

#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csCore/eigenUtils/RotationMatrix_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SO3.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

template <class Derived>
SO3<typename Derived::Scalar> so3Exp(Eigen::MatrixBase<Derived> const& delta) {
  using S = typename Derived::Scalar;
  CSLIE_ASSERT_DEBUG(
      csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SO3Alg<S>()),
      "delta required to be 3x1");

  SO3Alg<S> const w2 = delta.cwiseAbs2();
  S const theta_sq = w2.sum();  // sum of the squares
  S const theta = std::sqrt(theta_sq);

  S a, b;
  if (theta < Constants<S>::epsilon) {
    a = static_cast<S>(1.0) -
        theta_sq * (static_cast<S>(1.0 / 6.0) +
                    theta_sq * static_cast<S>(1.0 / 120.0));
    b = static_cast<S>(0.5) -
        theta_sq * (static_cast<S>(1.0 / 24.0) -
                    theta_sq * static_cast<S>(1.0 / 720.0));
  } else {
    a = sin(theta) / theta;
    b = (static_cast<S>(1.0) - cos(theta)) / theta_sq;
  }

  SO3<S> R;

  R(0, 0) = 1 - b * (w2(1) + w2(2));
  R(1, 1) = 1 - b * (w2(0) + w2(2));
  R(2, 2) = 1 - b * (w2(0) + w2(1));

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

template <class Derived>
SO3J<typename Derived::Scalar> so3J(Eigen::MatrixBase<Derived> const& delta) {
  using S = typename Derived::Scalar;
  CSLIE_ASSERT_DEBUG(
      csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SO3Alg<S>()),
      "delta required to be 3x1");

  S const th = delta.norm();

  if (th > Constants<S>::epsilon) {
    SO3J<S> ret;

    Eigen::Matrix<S, 3, 1> const a = delta.array() / th;

    const S sinth = std::sin(th);
    const S costh = std::cos(th);

    S tmp1 = 1.0 / th;
    S tmp2 = sinth * tmp1;
    S tmp3 = 1 - sinth * tmp1;
    S tmp4 = 1 - costh;
    S tmp5 = a(0) * a(1) * tmp3;
    S tmp6 = a(0) * a(2) * tmp3;
    S tmp7 = a(1) * a(2) * tmp3;

    ret(0, 0) = a(0) * a(0) * tmp3 + tmp2;
    ret(0, 1) = tmp5 - a(2) * tmp4 * tmp1;
    ret(0, 2) = tmp6 + a(1) * tmp4 * tmp1;
    ret(1, 0) = tmp5 + a(2) * tmp4 * tmp1;
    ret(1, 1) = a(1) * a(1) * tmp3 + tmp2;
    ret(1, 2) = tmp7 - a(0) * tmp4 * tmp1;
    ret(2, 0) = tmp6 - a(1) * tmp4 * tmp1;
    ret(2, 1) = tmp7 + a(0) * tmp4 * tmp1;
    ret(2, 2) = a(2) * a(2) * tmp3 + tmp2;

    return ret;
  } else {
    SO3J<S> ret;

    ret = Eigen::Matrix<S, 3, 3>::Identity() + 0.5 * skewSymm(delta);

    return ret;
  }
  // TODO(sc) why I cannot use UNREACHABLE here? it does not compile...
  CSLIE_ASSERT_PARANOID(
      "csLie::so3J this is a unreachable code, how did you ended up here?");
}

namespace detail {
template <typename T>
int signum(T val) {
  return (static_cast<T>(0) < val) - (val < static_cast<T>(0));
}
}  // namespace detail

// TODO(sc) references for better impl
// https://pettni.github.io/smooth/index.html
// https://github.com/strasdat/Sophus
// https://github.com/farm-ng/farm-ng-core/tree/cygnet-dev (aka Sophus 2.0)

template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 1> so3Log(
    Eigen::MatrixBase<Derived> const& rot) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(
      csCore::eigenUtils::isSameSizeStaticOrDynamic(rot, SO3<S>()),
      "rot required to be 3x3");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidRotationMatrix3D(rot),
                        "rot is not a valid rotation matrix");

  SO3Alg<S> ret;
  ret(0) = S(0.5) * (rot(2, 1) - rot(1, 2));
  ret(1) = S(0.5) * (rot(0, 2) - rot(2, 0));
  ret(2) = S(0.5) * (rot(1, 0) - rot(0, 1));

  S const ct = static_cast<S>(0.5) * (rot.trace() - 1);
  S const st2 = ret.cwiseAbs2().sum();
  S const st = std::sqrt(st2);

  // TODO(sc) find better constants / better impl
  static constexpr S almostOne{0.999856};
  static constexpr S almostMinusOne{-0.999856};

  if (ct > almostOne) {
    // special case at 1
    const S f = 1 + st2 * (static_cast<S>(1.0 / 6.0) +
                           st2 * (static_cast<S>(3.0 / 40.0) +
                                  st2 * static_cast<S>(5.0 / 112.0)));
    ret *= f;
  } else if (ct > almostMinusOne) {
    // standard case
    const S theta = std::acos(ct);
    const S f = theta / st;
    ret *= f;
  } else {
    // special case at -1
    const S theta = Constants<S>::pi - std::asin(st);
    const S invB = (theta * theta) / (1 - ct);

    const S w00 = invB * (rot(0, 0) - ct);
    const S w11 = invB * (rot(1, 1) - ct);
    const S w22 = invB * (rot(2, 2) - ct);

    const S w01 = invB * S(0.5) * (rot(1, 0) + rot(0, 1));
    const S w02 = invB * S(0.5) * (rot(2, 0) + rot(0, 2));
    const S w12 = invB * S(0.5) * (rot(2, 1) + rot(1, 2));

    if (w00 > w11) {
      if (w00 > w22) {
        ret(0) = detail::signum(ret(0)) * std::sqrt(w00);
        const S inv_w0 = 1.0 / ret(0);
        ret(1) = w01 * inv_w0;
        ret(2) = w02 * inv_w0;
      } else {
        ret(2) = detail::signum(ret(2)) * std::sqrt(w22);
        const S inv_w2 = 1 / ret(2);
        ret(0) = w02 * inv_w2;
        ret(1) = w12 * inv_w2;
      }
    } else if (w11 > w22) {
      ret(1) = detail::signum(ret(1)) * std::sqrt(w11);
      const S inv_w1 = 1 / ret(1);
      ret(0) = w01 * inv_w1;
      ret(2) = w12 * inv_w1;
    } else {
      ret(2) = detail::signum(ret(2)) * std::sqrt(w22);
      const S inv_w2 = 1 / ret(2);
      ret(0) = w02 * inv_w2;
      ret(1) = w12 * inv_w2;
    }
  }
  return ret;
}

template <class Derived>
SO3Jinv<typename Derived::Scalar> so3Jinv(
    Eigen::MatrixBase<Derived> const& delta) {
  using S = typename Derived::Scalar;
  CSLIE_ASSERT_DEBUG(
      csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, SO3Alg<S>()),
      "delta required to be 3x1");

  S const th = delta.norm();

  if (th > Constants<S>::epsilon) {
    SO3Jinv<S> ret;

    Eigen::Matrix<S, 3, 1> const a = delta.array() / th;

    S const sinth = std::sin(th);
    S const costh = std::cos(th);

    S const tmp1 = costh + 1;
    S const tmp2 = 1 / sinth;
    S const tmp3 = tmp1 * tmp2 * th / 2;
    S const tmp4 = 1 - tmp1 * tmp2 * th / 2;
    S const tmp5 = a(0) * a(1) * tmp4;
    S const tmp6 = a(0) * a(2) * tmp4;
    S const tmp7 = a(1) * a(2) * tmp4;

    ret(0, 0) = a(0) * a(0) * tmp4 + tmp3;
    ret(0, 1) = tmp5 + a(2) * th / 2;
    ret(0, 2) = tmp6 - a(1) * th / 2;
    ret(1, 0) = tmp5 - a(2) * th / 2;
    ret(1, 1) = a(1) * a(1) * tmp4 + tmp3;
    ret(1, 2) = tmp7 + a(0) * th / 2;
    ret(2, 0) = tmp6 + a(1) * th / 2;
    ret(2, 1) = tmp7 - a(0) * th / 2;
    ret(2, 2) = a(2) * a(2) * tmp4 + tmp3;

    return ret;
  } else {
    SO3Jinv<S> ret;

    ret = Eigen::Matrix<S, 3, 3>::Identity() - 0.5 * skewSymm(delta);

    return ret;
  }

  // TODO(sc) why I cannot use UNREACHABLE here? it does not compile...
  CSLIE_ASSERT_PARANOID(
      "csLie::so3Jinv this is a unreachable code, how did you ended up here?");
}

// reference for RPY conversion
// https://eecs.qmul.ac.uk/~gslabaugh/publications/euler.pdf
template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 1> so3ToRPY(
    Eigen::MatrixBase<Derived> const& rot) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(
      csCore::eigenUtils::isSameSizeStaticOrDynamic(rot, SO3<S>()),
      "rot required to be 3x3");
  CSLIE_ASSERT_PARANOID(csCore::eigenUtils::isValidRotationMatrix3D(rot),
                        "rot is not a valid rotation matrix");

  const S sin_pitch{rot(2, 0)};
  if (std::abs(std::abs(sin_pitch) - static_cast<S>(1.0)) >
      Constants<S>::epsilon) {
    const S pitch{-std::asin(rot(2, 0))};
    const S cos_pitch{std::cos(pitch)};  // can we extract it from the matrix?

    const S roll{std::atan2(rot(2, 1) / cos_pitch, rot(2, 2) / cos_pitch)};
    const S yaw{std::atan2(rot(1, 0) / cos_pitch, rot(0, 0) / cos_pitch)};
    return Eigen::Matrix<S, 3, 1>(roll, pitch, yaw);
  } else {
    if (sin_pitch < 0) {  // ~= -1
      const S pitch{Constants<S>::pi / static_cast<S>(2.0)};
      const S yaw = 0;
      const S roll = std::atan2(rot(0, 1), rot(0, 2));
      return Eigen::Matrix<S, 3, 1>(roll, pitch, yaw);
    } else {  // ~= 1
      const S pitch{-Constants<S>::pi / static_cast<S>(2.0)};
      const S yaw = 0;
      const S roll = std::atan2(-rot(0, 1), -rot(0, 2));
      return Eigen::Matrix<S, 3, 1>(roll, pitch, yaw);
    }
  }

  // TODO(sc) why I cannot use UNREACHABLE here? it does not compile...
  CSLIE_ASSERT_PARANOID(
      "csLie::so3ToRPY this is a unreachable code, how did you ended up here?");
}

template <class Derived>
SO3<typename Derived::Scalar> rpyToSo3(Eigen::MatrixBase<Derived> const& rpy) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(
                         rpy, Eigen::Matrix<S, 3, 1>()),
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

template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymm(
    Eigen::MatrixBase<Derived> const& v) {
  using S = typename Derived::Scalar;

  CSLIE_ASSERT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(
                         v, Eigen::Matrix<S, 3, 1>()),
                     "v required to be 3x1");

  Eigen::Matrix<S, 3, 3> R;

  R << 0, -v(2, 0), v(1, 0), v(2, 0), 0, -v(0, 0), -v(1, 0), v(0, 0), 0;

  return R;
}

}  // namespace csLie
