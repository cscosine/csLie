#pragma once

#include <Eigen/Core>
#include <csLie/Types.hpp>
#include <type_traits>

namespace csLie {

/// @brief compute the exp of a se2 algebra to a SE2
/// @tparam Derived template to support generic matrix
/// @param delta a 3d vector representing the se2 algebra vector
/// @return the rototranslation matrix in SE2 from the algebra vector in se2
/// @throw if delta is not a 3x1 matrix
template <class Derived>
SE2<typename Derived::Scalar> se2Exp(Eigen::MatrixBase<Derived> const& delta);

/*
TODO(sc) implement log in SE2
/// @brief check if a matrix is a valid rotation matrix, i.e., det(R) = 1, R*R^T
/// = I and R^T*R = I, withing a tolerance
/// @tparam Derived template to support generic matrix
/// @param rot the rotation matrix
/// @throw if rotation part of the transformation is not a 2x2 matrix (DEBUG
only), if rotation part of the transformation is  not a rotation
/// matrix within default precision (PARANOID)
template <class Derived>
typename Derived::Scalar se2Log(Eigen::MatrixBase<Derived> const& trans);
*/

/// @brief create a SE2 transformation from translation (x y) and angle (Y)
/// @tparam Derived template to support generic matrix
/// @param xyY a 3d vector representing the xy translation and the Y angle (YAW,
/// is around Z axis)
/// @return the rototranslation matrix in SE2
/// @throw if xyY is not a 3x1 matrix
template <class Derived>
SE2<typename Derived::Scalar> xyYTose2(Eigen::MatrixBase<Derived> const& xyY);

/// @brief extract translation (x y) and angle (Y) from a se2
/// @tparam Derived template to support generic matrix
/// @param transf a 3x3 matrix representing a isometry
/// @return the xyY vector, with a two element translation foloowed by the angle
/// @throw if transf is not a 3x3 valid isometry transformation
template <class Derived>
Eigen::Vector3<typename Derived::Scalar> se2ToxyY(
    Eigen::MatrixBase<Derived> const& transf);

/// @brief extract translation (x y) and angle (Y) from a se2
/// @tparam S the scalar type
/// @param transf a SE2, i.e., a Eigen::Transform<2, Isometry>
/// @return the xyY vector, with a two element translation foloowed by the angle
/// @throw if transf is not a 3x3 valid isometry transformation
template <class S>
Eigen::Vector3<S> se2ToxyY(SE2<S> const& transf);

}  // namespace csLie
