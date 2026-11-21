#pragma once

#include <Eigen/Core>
#include <csLie/Types.hpp>

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

/// @brief compute the RPY (roll pitch yaw) angles of a SO3 rotation
/// RPY are expressed as Euler angles so that R = rotZ(yaw) * rotY(pitch) *
/// rotY(roll)
/// Note: there are always two valid RPY for a rotation matrix,
/// this function returns the solution of so3ToRPYPair with pitch extracted as
/// asin(R(2,0)). In case a gymbal lock is detected, the solution with yaw=0.
/// is returned.
/// @tparam Derived template to support generic matrix
/// @param rot the rotation matrix
/// @throw if rot is not a 3x3 matrix (DEBUG only), if rot is not a rotation
/// matrix within default precision (PARANOID)
template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 1> so3ToRPY(Eigen::MatrixBase<Derived> const &rot);

/// @brief compute a SO3 rotation matrix from RPY (roll pitch yaw) angles
/// RPY are expressed as Euler angles so that R = rotZ(yaw) * rotY(pitch) *
/// rotY(roll)
/// @tparam Derived template to support generic matrix
/// @param rpy the roll pitch yaw angles
/// @throw if rpy is not a 3x1 matrix
template <class Derived> SO3<typename Derived::Scalar> rpyToSo3(Eigen::MatrixBase<Derived> const &rpy);

/// @brief return the skew symmetric matrix of a 3d vector
/// @tparam Derived template to support generic matrix
/// @param v the vector
/// @throw if v is not a 3x1 matrix
template <class Derived> Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymm(Eigen::MatrixBase<Derived> const &v);

} // namespace csLie
