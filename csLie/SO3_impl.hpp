#pragma once

#include <csCore/eigenUtils/MatrixShape_impl.hpp>
#include <csCore/eigenUtils/RotationMatrix_impl.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SO3.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

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

}  // namespace csLie
