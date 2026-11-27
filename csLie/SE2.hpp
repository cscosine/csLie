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
template <class Derived> SE2<typename Derived::Scalar> se2Exp(Eigen::MatrixBase<Derived> const &delta);

/// @brief compute the log of a SE2 rototranslation to a se2 algebra vector
/// @tparam Derived template to support generic matrix
/// @param transf a 3x3 matrix representing an isometry
/// @return the se2 algebra vector (x, y, theta), theta in (-pi, pi]: se2Exp of it gives transf
/// @throw if transf is not a 3x3 matrix (DEBUG only), if transf is not a valid isometry within default precision
/// (PARANOID)
template <class Derived> SE2Alg<typename Derived::Scalar> se2Log(Eigen::MatrixBase<Derived> const &transf);

/// @brief compute the log of a SE2 rototranslation to a se2 algebra vector
/// @tparam S the scalar type
/// @param transf a SE2, i.e., a Eigen::Transform<2, Isometry>
/// @return the se2 algebra vector (x, y, theta), theta in (-pi, pi]: se2Exp of it gives transf
/// @throw if transf is not a valid isometry within default precision (PARANOID)
template <class S> SE2Alg<S> se2Log(SE2<S> const &transf);

/// @brief create a SE2 rototranslation from a position (x, y) and a yaw angle in radians (around Z, counterclockwise):
/// R = so2Exp(yaw), t = (x, y). Not the se2 algebra: se2Exp((x, y, theta)) has a translation V (x, y) (an arc), equal
/// only at theta = 0
/// @tparam Derived template to support generic matrix
/// @param xyYaw a 3d vector (x, y, yaw), yaw in radians
/// @return the rototranslation in SE2
/// @throw if xyYaw is not a 3x1 matrix (DEBUG only)
template <class Derived> SE2<typename Derived::Scalar> xyYawToSe2(Eigen::MatrixBase<Derived> const &xyYaw);

/// @brief create a SE2 rototranslation from a position (x, y) and a yaw angle in radians, as xyYawToSe2(Vector3)
/// @tparam Derived template to support generic matrix
/// @param xy a 2d vector, the position
/// @param yawRad the yaw angle in radians
/// @return the rototranslation in SE2
/// @throw if xy is not a 2x1 matrix (DEBUG only)
template <class Derived>
SE2<typename Derived::Scalar> xyYawToSe2(Eigen::MatrixBase<Derived> const &xy, typename Derived::Scalar yawRad);

/// @brief create a SE2 rototranslation from a position (x, y) and a yaw angle in radians, as xyYawToSe2(Vector3)
/// @tparam S the scalar type
/// @param x, y the position
/// @param yawRad the yaw angle in radians
/// @return the rototranslation in SE2
template <class S> std::enable_if_t<std::is_floating_point<S>::value, SE2<S>> xyYawToSe2(S x, S y, S yawRad);

/// @brief the position (x, y) and the yaw angle in radians, in (-pi, pi], of a SE2 rototranslation: the inverse of
/// xyYawToSe2 (not the se2 log)
/// @tparam Derived template to support generic matrix
/// @param transf a 3x3 matrix representing an isometry
/// @return the vector (x, y, yaw)
/// @throw if transf is not a 3x3 matrix (DEBUG only), if transf is not a valid isometry (PARANOID)
template <class Derived> Eigen::Vector3<typename Derived::Scalar> se2ToXyYaw(Eigen::MatrixBase<Derived> const &transf);

/// @brief the position (x, y) and the yaw angle in radians, in (-pi, pi], of a SE2 rototranslation, as
/// se2ToXyYaw(matrix)
/// @tparam S the scalar type
/// @param transf a SE2, i.e., a Eigen::Transform<2, Isometry>
/// @return the vector (x, y, yaw)
/// @throw if transf is not a valid isometry (PARANOID)
template <class S> Eigen::Vector3<S> se2ToXyYaw(SE2<S> const &transf);

} // namespace csLie
