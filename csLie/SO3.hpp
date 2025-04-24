#pragma once

#include <Eigen/Core>
#include <csLie/Types.hpp>

namespace csLie {

/// @brief compute the RPY (roll pitch yaw) angles of a SO3 rotation
/// RPY are expressed as Euler angles so that R = rotZ(yaw) * rotY(pitch) *
/// rotY(roll)
/// Note: there are always two valid RPY for a rotation matrix,
/// this function returns the solution of so3ToRPYPair with pitch extracted as
/// asin(R(2,0)). In case a gymbal lock is detected, the solution with roll=0.
/// is returned.
/// @tparam Derived template to support generic matrix
/// @param rot the rotation matrix
/// @throw if rot is not a 3x3 matrix (DEBUG only), if rot is not a rotation
/// matrix within default precision (PARANOID)
template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 1> so3ToRPY(
    Eigen::MatrixBase<Derived> const& rot);

/// @brief compute a SO3 rotation matrxi from RPY (roll pitch yaw) angles
/// RPY are expressed as Euler angles so that R = rotZ(yaw) * rotY(pitch) *
/// rotY(roll)
/// @tparam Derived template to support generic matrix
/// @param rpy the roll pitch yaw angles
/// @throw if rpy is not a 3x1 matrix
template <class Derived>
SO3<typename Derived::Scalar> rpyToSo3(Eigen::MatrixBase<Derived> const& rpy);

}  // namespace csLie
