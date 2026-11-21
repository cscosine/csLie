#pragma once

#include <Eigen/Core>
#include <csLie/Types.hpp>

namespace csLie {

/// @brief compute the exp of a s2 algebra to a S2
/// @tparam Derived template to support generic matrix
/// @param delta a 2d vector representing the s2 algebra vector
/// @return the unit vector in S2 from the algebra vector in s2
template <class Derived> S2<typename Derived::Scalar> s2Exp(Eigen::MatrixBase<Derived> const &delta);

/// @brief transform a S2 unit vector into a SO3 rotation matrix with x aligned
/// to the unit vector
/// @tparam Derived template to support generic matrix
/// @param v the unit vector
/// @return a rotation matrix
template <class Derived> SO3<typename Derived::Scalar> s2toZAlignedSO3(Eigen::MatrixBase<Derived> const &v);

/// @brief perform the "oplus" addition of a S2 element with the exp of a s2
/// algebra vector given
/// - a unitary 3D vector v (in S2)
/// - a s2 algebra vector delta (2D vector)
/// return the oplus(v, delta)
/// equivalent to s2toZAlignedSO3(v) * s2Exp(delta)
/// for readability rewritten as
/// R(V) * exp(delta)
/// note: repeating oplus in loops might introduce numerical errors, renormalize
/// the vector when needed.
///
/// @tparam Derived1 template to support generic matrix
/// @tparam Derived2 template to support generic matrix
/// @param v a 3d unit vector representing the s2
/// @param delta a 2d vector representing the s2 algebra vector
/// @return the result of the oplus, namely s2toZAlignedSO3(V) * s2Exp(delta)
template <class Derived1, class Derived2>
S2<typename Derived1::Scalar> s2oplus(Eigen::MatrixBase<Derived1> const &v, Eigen::MatrixBase<Derived2> const &delta);

/// @brief compute the jacobian of the oplus operation when delta is zero
///
/// the typical operation `oplus` given v in S2 and delta as S2Alg is defined as
/// v_p = s2toZAlignedSO3(v) * s2Exp(delta)
/// for readability rewritten as
/// v_p = R(V) * exp(delta)
/// in an optimization framework, is often need to evaluate the jacobian wrt the
/// delta evaluated at zero
///  d v_p  |
/// ------- |
/// d delta |delta=0
/// this is equal to the last two column of R(v), which are returned by this
/// function
/// ///
/// @tparam Derived the matrix derived type
/// @param v the S2 element
/// @return the jacobians of the oplus(v, delta) evaluated in delta = 0
template <class Derived> s2oplusJac<typename Derived::Scalar> s2_d_oplus_d_delta_0(Eigen::MatrixBase<Derived> const &v);

} // namespace csLie
