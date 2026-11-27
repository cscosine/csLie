#pragma once

#include <Eigen/Core>
#include <csLie/Constants.hpp>
#include <csLie/Types.hpp>
#include <optional>

namespace csLie {

/// @brief compute the exp of a so3 algebra to a SO3
/// @tparam Derived template to support generic matrix
/// @param delta a 3d vector representing the so3 algebra vector
/// @return the rotation matrix in SO3 from the algebra vector in so3
/// @throw if delta is not a 3x1 matrix
template <class Derived> SO3<typename Derived::Scalar> so3Exp(Eigen::MatrixBase<Derived> const &delta);

/// @brief compute the J, i.e. the left jacobian of a so3 algebra
/// @tparam Derived template to support generic matrix
/// @param delta a 3d vector representing the so3 algebra vector
/// @return the J matrix from the algebra vector in so3
/// @throw if delta is not a 3x1 matrix
template <class Derived> SO3J<typename Derived::Scalar> so3J(Eigen::MatrixBase<Derived> const &delta);

/// @brief compute the log of a SO3 rotation to a so3 vector
/// @tparam Derived template to support generic matrix
/// @param rot the rotation matrix
/// @throw if rot is not a 3x3 matrix (DEBUG only), if rot is not a rotation
/// matrix within default precision (PARANOID)
template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 1> so3Log(Eigen::MatrixBase<Derived> const &rot);

/// @brief compute the J^-1, i.e. the inverse left jacobian of a so3 algebra
/// @tparam Derived template to support generic matrix
/// @param delta a 3d vector representing the so3 algebra vector
/// @return the J^-1 matrix from the algebra vector in so3
/// @throw if delta is not a 3x1 matrix
template <class Derived> SO3Jinv<typename Derived::Scalar> so3Jinv(Eigen::MatrixBase<Derived> const &delta);

/// @brief where so3ToRpy puts the free angle at gimbal lock (pitch = +-pi / 2, where only yaw - roll at pitch = pi / 2,
/// yaw + roll at pitch = -pi / 2, is determined)
enum class GimbalPolicy {
  kRollZero, ///< roll = 0, the free angle in yaw
  kYawZero,  ///< yaw = 0, the free angle in roll
};

/// @brief the default gimbal tolerance of so3ToRpyChecked: sqrt(epsilon). Roll and yaw, each alone, are accurate to
/// about epsilon / cos(pitch): at cos(pitch) = sqrt(epsilon), half their digits
template <class S> S rpyCheckedGimbalTolerance();

/// @brief compute the RPY (roll, pitch, yaw) angles of a SO3 rotation, in radians, R = rotZ(yaw) rotY(pitch)
/// rotX(roll): pitch in [-pi / 2, pi / 2], roll and yaw in (-pi, pi]. rpyToSo3 of the result gives rot back to a few
/// epsilon at every pitch. Roll and yaw, each alone, are accurate to about epsilon / cos(pitch); at gimbal lock,
/// cos(pitch) <= gimbalTolerance, only yaw - roll (pitch = pi / 2) or yaw + roll (pitch = -pi / 2) is determined: the
/// angle the policy names is 0 and the other carries the rotation. That drops terms of size cos(pitch) from the
/// rotation: a gimbalTolerance larger than the default epsilon costs a round trip error up to it
/// @tparam P where the free angle goes at gimbal lock (default: roll = 0)
/// @tparam Derived template to support generic matrix
/// @param rot the rotation matrix
/// @param gimbalTolerance the cos(pitch) at and below which the rotation is treated as at gimbal lock
/// @throw if rot is not a 3x3 matrix (DEBUG only), if rot is not a rotation matrix within default precision (PARANOID)
template <GimbalPolicy P = GimbalPolicy::kRollZero, class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 1>
so3ToRpy(Eigen::MatrixBase<Derived> const &rot,
         typename Derived::Scalar gimbalTolerance = Constants<typename Derived::Scalar>::epsilon);

/// @brief as so3ToRpy, empty when roll and yaw, each alone, are not accurate: cos(pitch) <= gimbalTolerance (default
/// sqrt(epsilon): half their digits)
/// @tparam Derived template to support generic matrix
/// @param rot the rotation matrix
/// @param gimbalTolerance the cos(pitch) at and below which the result is empty
/// @throw if rot is not a 3x3 matrix (DEBUG only), if rot is not a rotation matrix within default precision (PARANOID)
template <class Derived>
std::optional<Eigen::Matrix<typename Derived::Scalar, 3, 1>>
so3ToRpyChecked(Eigen::MatrixBase<Derived> const &rot,
                typename Derived::Scalar gimbalTolerance = rpyCheckedGimbalTolerance<typename Derived::Scalar>());

/// @brief compute a SO3 rotation matrix from RPY (roll, pitch, yaw) angles in radians, R = rotZ(yaw) rotY(pitch)
/// rotX(roll)
/// @tparam Derived template to support generic matrix
/// @param rpy the roll, pitch and yaw angles
/// @throw if rpy is not a 3x1 matrix
template <class Derived> SO3<typename Derived::Scalar> rpyToSo3(Eigen::MatrixBase<Derived> const &rpy);

/// @brief return the skew symmetric matrix of a 3d vector
/// @tparam Derived template to support generic matrix
/// @param v the vector
/// @throw if v is not a 3x1 matrix
template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymm(Eigen::MatrixBase<Derived> const &v);

} // namespace csLie
