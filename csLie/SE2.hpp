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

/// @brief create a SE2 transformation from translation (x y) and angle (Y)
/// @tparam Derived template to support generic matrix
/// @param xyY a 3d vector representing the xy translation and the Y angle (YAW,
/// is around Z axis)
/// @return the rototranslation matrix in SE2
/// @throw if xyY is not a 3x1 matrix
template <class Derived> SE2<typename Derived::Scalar> xyYToSe2(Eigen::MatrixBase<Derived> const &xyY);

/// @brief extract translation (x y) and angle (Y) from a se2
/// @tparam Derived template to support generic matrix
/// @param transf a 3x3 matrix representing a isometry
/// @return the xyY vector, with a two element translation followed by the angle
/// @throw if transf is not a 3x3 valid isometry transformation
template <class Derived> Eigen::Vector3<typename Derived::Scalar> se2ToxyY(Eigen::MatrixBase<Derived> const &transf);

/// @brief extract translation (x y) and angle (Y) from a se2
/// @tparam S the scalar type
/// @param transf a SE2, i.e., a Eigen::Transform<2, Isometry>
/// @return the xyY vector, with a two element translation followed by the angle
/// @throw if transf is not a 3x3 valid isometry transformation
template <class S> Eigen::Vector3<S> se2ToxyY(SE2<S> const &transf);

} // namespace csLie
