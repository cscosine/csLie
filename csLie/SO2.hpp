#pragma once

#include <Eigen/Core>

#include <csCore/core/BasicTypes.hpp>
#include <csLie/Types.hpp>
#include <csLie/Constants.hpp>

#include <type_traits>

namespace csLie {

  /// @brief compute the exp of a so2 algebra to a SO2, version with return type
  /// @tparam S scalar type
  /// @param delta a 1d vector representing the so2 algebra vector
  /// @return the rotation matrix in SO2 from the algebra vector in so2
  /// @throw if delta is not a 1x1 matrix
  template <class Derived> 
  SO2<typename Derived::Scalar> so2Exp(Eigen::MatrixBase<Derived> const& delta);

  /// @brief compute the exp of a so2 algebra to a SO2, version with scalar input and return type
  /// @tparam S scalar type
  /// @param delta a scalar vector representing the so2 algebra vector
  /// @return the rotation matrix in SO2 from the algebra vector in so2
  template <class S, class RetType = SO2<S>>
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta);

  /// @brief check if a matrix is a valid rotation matrix, i.e., det(R) = 1, R*R^T = I and R^T*R = I, withing a tolerance
  /// @param rot the rotation matrix
  /// @throw if rot is not a 2x2 matrix (DEBUG only), if rot is not a rotation matrix within default precision (PARANOID)
  template <class Derived>
  typename Derived::Scalar so2Log(Eigen::MatrixBase<Derived> const& rot);

}