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

}  // namespace csLie
