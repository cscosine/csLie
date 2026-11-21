#pragma once

#include <Eigen/Core>
#include <csLie/Types.hpp>

namespace csLie {

/// @brief normalize a 3d vector and compute its jacobian
/// @tparam Derived template to support generic matrix
/// @param v a 3d vector
/// @return the unit vector vo and its jacobian d vo / d v
/// @throw if v is not a 3x1 matrix (DEBUG)
/// @throw if v has zero norm and cannot be normalized (PARANOID)
template <class Derived>
UnitVector3dWithJacobians<typename Derived::Scalar> normalize3dVectorWithJacobians(Eigen::MatrixBase<Derived> const &v);

// TODO(sc) implement a version returning a optional in case the normalization
// cannot be performed

} // namespace csLie
