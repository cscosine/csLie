#pragma once

#include <Eigen/Core>
#include <csLie/Types.hpp>

namespace csLie {

/// @brief compute the exp of a se3 algebra to a SE3
/// @tparam Derived template to support generic matrix
/// @param delta a 6d vector representing the se3 algebra vector
/// @return the matrix in SE3 from the algebra vector in se3
/// @throw if delta is not a 6x1 matrix
template <class Derived>
SE3<typename Derived::Scalar> se3Exp(Eigen::MatrixBase<Derived> const& delta);

/// @brief compute the log of a SE3 rotation to a se3 vector
/// @tparam Derived template to support generic matrix
/// @param transf the transformation matrix
/// @throw if rot is not a 3x3 matrix (DEBUG only), if rot is not a isometry
/// matrix within default precision (PARANOID)
template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 6, 1> se3Log(
    Eigen::MatrixBase<Derived> const& transf);

}  // namespace csLie
