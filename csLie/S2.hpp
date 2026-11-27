#pragma once

#include <Eigen/Core>
#include <csLie/Types.hpp>

namespace csLie {

/// @brief the exp of S2 at e_x = (1, 0, 0): (cos(theta), sin(theta) / theta delta), theta = |delta|, the point at the
/// angle theta from e_x in the direction delta of its tangent plane (the y, z plane)
/// @tparam Derived template to support generic matrix
/// @param delta a 2d vector, the S2 tangent vector at e_x
/// @return the unit vector
/// @throw if delta is not a 2x1 matrix (DEBUG only)
template <class Derived> S2<typename Derived::Scalar> s2Exp(Eigen::MatrixBase<Derived> const &delta);

/// @brief the frame of a unit vector, the rotation [v, b1, b2] that takes e_x to v; (b1, b2) is the basis of the
/// tangent plane at v of s2Oplus and s2OplusInverse. b2 = (a x v) / |a x v|, b1 = b2 x v with the pivot axis a = e_x
/// while v_x^2 <= 1/2, else e_z: two charts, both well conditioned (|a x v| >= 1/sqrt(2)), no trigonometry. The frame
/// jumps between them (every frame of S2 is discontinuous somewhere)
/// @tparam Derived template to support generic matrix
/// @param v the unit vector
/// @return the rotation matrix [v, b1, b2]
/// @throw if v is not a 3x1 matrix (DEBUG only), if v is not a unit vector within default precision (PARANOID)
template <class Derived> SO3<typename Derived::Scalar> s2Frame(Eigen::MatrixBase<Derived> const &v);

/// @brief the "oplus" of S2 (g2o's word): v moved by delta in its tangent plane, s2Frame(v) s2Exp(delta). S2 is not a
/// group: one oplus, a chart at v through s2Frame(v), no left or right. Repeated in a loop, renormalize the result
/// when needed
/// @tparam Derived1 template to support generic matrix
/// @tparam Derived2 template to support generic matrix
/// @param v the unit vector
/// @param delta a 2d vector, the S2 tangent vector in the basis (b1, b2) of s2Frame(v)
/// @return the unit vector s2Frame(v) s2Exp(delta)
/// @throw if v is not a 3x1 or delta not a 2x1 matrix (DEBUG only), if v is not a unit vector (PARANOID)
template <class Derived1, class Derived2>
S2<typename Derived1::Scalar> s2Oplus(Eigen::MatrixBase<Derived1> const &v, Eigen::MatrixBase<Derived2> const &delta);

/// @brief the inverse of s2Oplus at v: the delta, in the basis (b1, b2) of s2Frame(v), with s2Oplus(v, delta) = w and
/// |delta| in [0, pi]. At w = -v every direction gives w: the one of b1, delta = (pi, 0); near -v the direction of
/// delta jumps (the cut locus), s2Oplus(v, delta) = w holds everywhere
/// @tparam Derived1 template to support generic matrix
/// @tparam Derived2 template to support generic matrix
/// @param v the unit vector whose chart is used
/// @param w the unit vector to reach
/// @return the S2 tangent vector delta
/// @throw if v or w is not a 3x1 matrix (DEBUG only), if v is not a unit vector (PARANOID)
template <class Derived1, class Derived2>
S2Tangent<typename Derived1::Scalar> s2OplusInverse(Eigen::MatrixBase<Derived1> const &v,
                                                    Eigen::MatrixBase<Derived2> const &w);

/// @brief the jacobian of s2Oplus(v, delta) in delta at delta = 0: the basis (b1, b2), the last two columns of
/// s2Frame(v)
/// @tparam Derived template to support generic matrix
/// @param v the unit vector
/// @return the 3x2 jacobian
/// @throw if v is not a 3x1 matrix (DEBUG only), if v is not a unit vector (PARANOID)
template <class Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 2> s2OplusJacobianAtZero(Eigen::MatrixBase<Derived> const &v);

} // namespace csLie
