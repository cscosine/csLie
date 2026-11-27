#pragma once

#include <Eigen/Core>
#include <csLie/Constants.hpp>
#include <csLie/Types.hpp>
#include <optional>

namespace csLie {

/// @brief compute the exp of a so3 tangent vector to a SO3
/// @tparam Derived template to support generic matrix
/// @param delta a 3d vector representing the so3 tangent vector
/// @return the rotation matrix in SO3 from the so3 tangent vector
/// @throw if delta is not a 3x1 matrix
template <class Derived> SO3<typename Derived::Scalar> so3Exp(Eigen::MatrixBase<Derived> const &delta);

/// @brief the left jacobian of SO3, Jl(delta): exp(delta + d) = exp(Jl(delta) d) exp(delta) to first order in d
/// @tparam Derived template to support generic matrix
/// @param delta a 3d vector representing the so3 tangent vector
/// @return the 3x3 left jacobian
/// @throw if delta is not a 3x1 matrix
template <class Derived>
Eigen::Matrix3<typename Derived::Scalar> so3LeftJacobian(Eigen::MatrixBase<Derived> const &delta);

/// @brief compute the log of a SO3 rotation to a so3 vector
/// @tparam Derived template to support generic matrix
/// @param rot the rotation matrix
/// @throw if rot is not a 3x3 matrix (DEBUG only), if rot is not a rotation
/// matrix within default precision (PARANOID)
template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 1> so3Log(Eigen::MatrixBase<Derived> const &rot);

/// @brief the inverse of the left jacobian of SO3, Jl(delta)^-1 (singular at |delta| = 2 pi)
/// @tparam Derived template to support generic matrix
/// @param delta a 3d vector representing the so3 tangent vector
/// @return the 3x3 inverse left jacobian
/// @throw if delta is not a 3x1 matrix
template <class Derived>
Eigen::Matrix3<typename Derived::Scalar> so3LeftJacobianInverse(Eigen::MatrixBase<Derived> const &delta);

/// @brief where so3ToRollPitchYaw puts the free angle at gimbal lock (pitch = +-pi / 2, where only yaw - roll at pitch
/// = pi / 2, yaw + roll at pitch = -pi / 2, is determined)
enum class GimbalPolicy {
  kRollZero, ///< roll = 0, the free angle in yaw
  kYawZero,  ///< yaw = 0, the free angle in roll
};

/// @brief the default gimbal tolerance of so3ToRollPitchYawChecked: sqrt(epsilon). Roll and yaw, each alone, are
/// accurate to about epsilon / cos(pitch): at cos(pitch) = sqrt(epsilon), half their digits
template <class S> S rollPitchYawCheckedGimbalTolerance();

/// @brief compute the RPY (roll, pitch, yaw) angles of a SO3 rotation, in radians, R = rotZ(yaw) rotY(pitch)
/// rotX(roll): pitch in [-pi / 2, pi / 2], roll and yaw in (-pi, pi]. rollPitchYawToSo3 of the result gives rot back to
/// a few epsilon at every pitch. Roll and yaw, each alone, are accurate to about epsilon / cos(pitch); at gimbal lock,
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
so3ToRollPitchYaw(Eigen::MatrixBase<Derived> const &rot,
                  typename Derived::Scalar gimbalTolerance = Constants<typename Derived::Scalar>::epsilon);

/// @brief as so3ToRollPitchYaw, empty when roll and yaw, each alone, are not accurate: cos(pitch) <= gimbalTolerance
/// (default sqrt(epsilon): half their digits)
/// @tparam Derived template to support generic matrix
/// @param rot the rotation matrix
/// @param gimbalTolerance the cos(pitch) at and below which the result is empty
/// @throw if rot is not a 3x3 matrix (DEBUG only), if rot is not a rotation matrix within default precision (PARANOID)
template <class Derived>
std::optional<Eigen::Matrix<typename Derived::Scalar, 3, 1>> so3ToRollPitchYawChecked(
    Eigen::MatrixBase<Derived> const &rot,
    typename Derived::Scalar gimbalTolerance = rollPitchYawCheckedGimbalTolerance<typename Derived::Scalar>());

/// @brief compute a SO3 rotation matrix from RPY (roll, pitch, yaw) angles in radians, R = rotZ(yaw) rotY(pitch)
/// rotX(roll)
/// @tparam Derived template to support generic matrix
/// @param rpy the roll, pitch and yaw angles
/// @throw if rpy is not a 3x1 matrix
template <class Derived> SO3<typename Derived::Scalar> rollPitchYawToSo3(Eigen::MatrixBase<Derived> const &rpy);

/// @brief the right jacobian of SO3, Jr(delta) = Jl(-delta) = Jl(delta)^T: exp(delta + d) = exp(delta) exp(Jr(delta) d)
/// to first order in d
/// @tparam Derived template to support generic matrix
/// @param delta a 3d vector, the so3 tangent vector
/// @return the 3x3 right jacobian
/// @throw if delta is not a 3x1 matrix (DEBUG only)
template <class Derived>
Eigen::Matrix3<typename Derived::Scalar> so3RightJacobian(Eigen::MatrixBase<Derived> const &delta);

/// @brief the inverse of the right jacobian of SO3, Jr(delta)^-1 = Jl(delta)^-T (singular at |delta| = 2 pi)
/// @tparam Derived template to support generic matrix
/// @param delta a 3d vector, the so3 tangent vector
/// @return the 3x3 inverse right jacobian
/// @throw if delta is not a 3x1 matrix (DEBUG only)
template <class Derived>
Eigen::Matrix3<typename Derived::Scalar> so3RightJacobianInverse(Eigen::MatrixBase<Derived> const &delta);

/// @brief the rotation perturbed on the left, exp(delta) rot (delta in the world frame)
/// @tparam Derived1, Derived2 template to support generic matrix
/// @param rot the rotation matrix
/// @param delta a 3d vector, the so3 tangent vector
/// @throw if rot is not a 3x3 or delta not a 3x1 matrix (DEBUG only)
template <class Derived1, class Derived2>
SO3<typename Derived1::Scalar> so3OplusLeft(Eigen::MatrixBase<Derived1> const &rot,
                                            Eigen::MatrixBase<Derived2> const &delta);

/// @brief the rotation perturbed on the right, rot exp(delta) (delta in the body frame)
/// @tparam Derived1, Derived2 template to support generic matrix
/// @param rot the rotation matrix
/// @param delta a 3d vector, the so3 tangent vector
/// @throw if rot is not a 3x3 or delta not a 3x1 matrix (DEBUG only)
template <class Derived1, class Derived2>
SO3<typename Derived1::Scalar> so3OplusRight(Eigen::MatrixBase<Derived1> const &rot,
                                             Eigen::MatrixBase<Derived2> const &delta);

/// @brief the skew symmetric matrix of a 3d vector (so3 hat): hat(v) w = v x w
/// @tparam Derived template to support generic matrix
/// @param v the vector
/// @throw if v is not a 3x1 matrix (DEBUG only)
template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 3> hat(Eigen::MatrixBase<Derived> const &v);

/// @brief the 3d vector of a skew symmetric matrix (so3 vee), the inverse of hat: (W(2, 1), W(0, 2), W(1, 0))
/// @tparam Derived template to support generic matrix
/// @param W the skew symmetric matrix (not checked)
/// @throw if W is not a 3x3 matrix (DEBUG only)
template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 1> vee(Eigen::MatrixBase<Derived> const &W);

} // namespace csLie
