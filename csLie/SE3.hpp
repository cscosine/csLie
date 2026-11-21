#pragma once

#include <Eigen/Core>
#include <csLie/Types.hpp>

namespace csLie {

/// @brief compute the exp of a se3 algebra to a SE3
/// @tparam Derived template to support generic matrix
/// @param delta a 6d vector representing the se3 algebra vector
/// @return the matrix in SE3 from the algebra vector in se3
/// @throw if delta is not a 6x1 matrix
template <class Derived> SE3<typename Derived::Scalar> se3Exp(Eigen::MatrixBase<Derived> const &delta);

/// @brief compute the log of a SE3 transformation to a se3 vector
/// @tparam Derived template to support generic matrix
/// @param transf the transformation matrix
/// @throw if transf is not a 4x4 matrix (DEBUG only), if transf is not a
/// isometry matrix within default precision (PARANOID)
template <class Derived> SE3Alg<typename Derived::Scalar> se3Log(Eigen::MatrixBase<Derived> const &transf);

/// @brief se3Log, but accepting a SE3<S> i.e., a Eigen::Transform<3, Isometry>
/// @tparam S the scalar type
/// @param transf a SE3, i.e., a Eigen::Transform<3, Isometry>
/// @throw if transf is not a 3x3 valid isometry transformation
template <class S> SE3Alg<S> se3Log(SE3<S> const &transf);

/// @brief compute the adjoint of a SE3 transformation
/// @tparam Derived template to support generic matrix
/// @param transf the transformation matrix
/// @throw if transf is not a 4x4 matrix (DEBUG only), if transf is not a
/// isometry matrix within default precision (PARANOID)
template <class Derived> SE3Adj<typename Derived::Scalar> se3Adjoint(Eigen::MatrixBase<Derived> const &transf);

/// @brief se3Adjoint, but accepting a SE3<S> i.e., a Eigen::Transform<3,
/// Isometry>
template <class S> SE3Adj<S> se3Adjoint(SE3<S> const &transf);

/// @brief SE3 d exp(eps)*T*p / d eps |_ eps=0
/// @tparam Derived template to support generic matrix
/// @param transformedPoint the 3x1 3d point transformed
/// @throw if transformedPoint is not a 3x1 matrix (DEBUG only)
template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 6>
se3JacobianLeftZeroPerturbationPointTransformation(Eigen::MatrixBase<Derived> const &transformedPoint);

} // namespace csLie
