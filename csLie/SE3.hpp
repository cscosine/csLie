#pragma once

#include <Eigen/Core>
#include <csLie/Types.hpp>

namespace csLie {

/// @brief compute the exp of a se3 tangent vector to a SE3
/// @tparam Derived template to support generic matrix
/// @param delta a 6d vector representing the se3 tangent vector
/// @return the matrix in SE3 from the se3 tangent vector
/// @throw if delta is not a 6x1 matrix
template <class Derived> SE3<typename Derived::Scalar> se3Exp(Eigen::MatrixBase<Derived> const &delta);

/// @brief compute the log of a SE3 transformation to a se3 vector
/// @tparam Derived template to support generic matrix
/// @param transf the transformation matrix
/// @throw if transf is not a 4x4 matrix (DEBUG only), if transf is not a
/// isometry matrix within default precision (PARANOID)
template <class Derived> SE3Tangent<typename Derived::Scalar> se3Log(Eigen::MatrixBase<Derived> const &transf);

/// @brief se3Log, but accepting a SE3<S> i.e., a Eigen::Transform<3, Isometry>
/// @tparam S the scalar type
/// @param transf a SE3, i.e., a Eigen::Transform<3, Isometry>
/// @throw if transf is not a 3x3 valid isometry transformation
template <class S> SE3Tangent<S> se3Log(SE3<S> const &transf);

/// @brief compute the adjoint of a SE3 transformation
/// @tparam Derived template to support generic matrix
/// @param transf the transformation matrix
/// @throw if transf is not a 4x4 matrix (DEBUG only), if transf is not a
/// isometry matrix within default precision (PARANOID)
template <class Derived> SE3Adjoint<typename Derived::Scalar> se3Adjoint(Eigen::MatrixBase<Derived> const &transf);

/// @brief se3Adjoint, but accepting a SE3<S> i.e., a Eigen::Transform<3,
/// Isometry>
template <class S> SE3Adjoint<S> se3Adjoint(SE3<S> const &transf);

/// @brief the jacobian of a transformed point under a left perturbation, d (exp(eps) T p) / d eps at eps = 0, eps =
/// (rho, omega) in the world frame: [I, -hat(T p)]
/// @tparam Derived template to support generic matrix
/// @param transformedPoint the 3x1 point already transformed, T p
/// @throw if transformedPoint is not a 3x1 matrix (DEBUG only)
template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 6> se3PointJacobianLeft(Eigen::MatrixBase<Derived> const &transformedPoint);

/// @brief the jacobian of a transformed point under a right perturbation, d (T exp(eps) p) / d eps at eps = 0, eps =
/// (rho, omega) in the body frame: R [I, -hat(p)]
/// @tparam Derived1, Derived2 template to support generic matrix
/// @param rot the rotation of T
/// @param point the 3x1 point, not transformed
/// @throw if rot is not a 3x3 or point not a 3x1 matrix (DEBUG only)
template <class Derived1, class Derived2>
Eigen::Matrix<typename Derived1::Scalar, 3, 6> se3PointJacobianRight(Eigen::MatrixBase<Derived1> const &rot,
                                                                     Eigen::MatrixBase<Derived2> const &point);

/// @brief the transformation perturbed on the left, exp(delta) T (delta = (rho, omega) in the world frame)
/// @tparam S the scalar type
/// @tparam Derived template to support generic matrix
/// @param transf the transformation
/// @param delta a 6d vector, the se3 tangent vector
/// @throw if delta is not a 6x1 matrix (DEBUG only)
template <class S, class Derived> SE3<S> se3OplusLeft(SE3<S> const &transf, Eigen::MatrixBase<Derived> const &delta);

/// @brief the transformation perturbed on the right, T exp(delta) (delta = (rho, omega) in the body frame)
/// @tparam S the scalar type
/// @tparam Derived template to support generic matrix
/// @param transf the transformation
/// @param delta a 6d vector, the se3 tangent vector
/// @throw if delta is not a 6x1 matrix (DEBUG only)
template <class S, class Derived> SE3<S> se3OplusRight(SE3<S> const &transf, Eigen::MatrixBase<Derived> const &delta);

} // namespace csLie
