#pragma once

#include <Eigen/Core>

#include <csCore/core/BasicTypes.hpp>
#include <csLie/Types.hpp>

#include <type_traits>

namespace csLie {

  // --- delta type matrix ---

  /// @brief compute the exp of a so2 algebra to a SO2, version with return type
  /// @tparam S scalar type
  /// @param delta a 1d vector representing the so2 algebra vector
  /// @return the rotation matrix in SO2 from the algebra vector in so2
  template <class Derived> 
  SO2<typename Derived::Scalar> so2Exp(Eigen::MatrixBase<Derived> const& delta);

  /// @brief compute the exp of a so2 algebra to a SO2, storing result on a output parameter
  /// @tparam S scalar type
  /// @param delta a 1d vector representing the so2 algebra vector
  /// @param res the rotation matrix in SO2 from the algebra vector in so2
  template <class DerivedIn, class DerivedOut> 
  void so2Exp(Eigen::MatrixBase<DerivedIn> const& delta, Eigen::MatrixBase<DerivedOut> & res);

  /// @brief compute the exp of a so2 algebra to a SO2, storing result on a output block
  /// @tparam S scalar type
  /// @param delta a 1d vector representing the so2 algebra vector
  /// @param res the rotation matrix in SO2 from the algebra vector in so2
  template <class DerivedIn, class DerivedOut, int BlockRows, int BlockCols, bool InnerPanel> 
  void so2Exp(Eigen::MatrixBase<DerivedIn> const& delta, Eigen::Block<DerivedOut, BlockRows, BlockCols, InnerPanel> res);
  
  // --- delta type scalar ---

  /// @brief compute the exp of a so2 algebra to a SO2, version with scalar input and return type
  /// @tparam S scalar type
  /// @param delta a scalar vector representing the so2 algebra vector
  /// @return the rotation matrix in SO2 from the algebra vector in so2
  template <class S, class RetType = SO2<S>>
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta);

  /// @brief compute the exp of a so2 algebra to a SO2, storing result on a output parameter
  /// @tparam S scalar type
  /// @param delta a scalar representing the so2 algebra vector
  /// @param res the rotation matrix in SO2 from the algebra vector in so2
  template <class S, class DerivedOut, class RetType = void> 
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta, Eigen::MatrixBase<DerivedOut>& res);

  /// @brief compute the exp of a so2 algebra to a SO2, storing result on a block
  /// @tparam S scalar type
  /// @param delta scalar representing the so2 algebra vector
  /// @param res the rotation matrix in SO2 from the algebra vector in so2
  template <class S, class DerivedOut, int BlockRows, int BlockCols, bool InnerPanel, class RetType = void> 
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta, Eigen::Block<DerivedOut, BlockRows, BlockCols, InnerPanel> res);

  // --- delta type block ---

  /// @brief compute the exp of a so2 algebra to a SO2, version with block input and return type
  /// @tparam S scalar type
  /// @param delta a 1x1 block representing the so2 algebra vector
  /// @return the rotation matrix in SO2 from the algebra vector in so2
  template <class DerivedIn, int BlockRows, int BlockCols, bool InnerPanel> 
  SO2<typename DerivedIn::Scalar> so2Exp(Eigen::Block<DerivedIn, BlockRows, BlockCols, InnerPanel> const& delta);

  /// @brief compute the exp of a so2 algebra to a SO2, storing result on a output parameter
  /// @tparam S scalar type
  /// @param delta a 1x1 block representing the so2 algebra vector
  /// @param res the rotation matrix in SO2 from the algebra vector in so2
  template <class DerivedIn, int BlockRows, int BlockCols, bool InnerPanel, class DerivedOut> 
  void so2Exp(Eigen::Block<DerivedIn, BlockRows, BlockCols, InnerPanel> const& delta, Eigen::MatrixBase<DerivedOut> & res);

  /// @brief compute the exp of a so2 algebra to a SO2, storing result on a block
  /// @tparam S scalar type
  /// @param delta a 1x1 block representing the so2 algebra vector
  /// @param res the rotation matrix in SO2 from the algebra vector in so2
  template <class DerivedIn, int BlockRowsIn, int BlockColsIn, bool InnerPanelIn, class DerivedOut, int BlockRowsOut, int BlockColsOut, bool InnerPanelOut> 
  void so2Exp(Eigen::Block<DerivedIn, BlockRowsIn, BlockColsIn, InnerPanelIn> const& delta, Eigen::Block<DerivedOut, BlockRowsOut, BlockColsOut, InnerPanelOut> res);

}