#pragma once

#include <Eigen/Core>

#include <csCore/BasicTypes.hpp>
#include <csLie/Types.hpp>

#include <type_traits>

namespace csLie {

  // TODO(sceriani) actualize
  /// @brief compute the exp of a so2 algebra to a SO2
  /// @tparam S scalar type
  /// @param delta a 1d vector representing the so2 algebra vector
  /// @return the rotation matrix in SO2 from the algebra vector in so2
  //template <class S> SO2<S> so2Exp(SO2Alg<S> const& delta);
  template <class Derived> SO2<typename Derived::Scalar> so2Exp(Eigen::MatrixBase<Derived> const& delta);

  template <class DerivedIn, class DerivedOut> void so2Exp(Eigen::MatrixBase<DerivedIn> const& delta, Eigen::MatrixBase<DerivedOut> & res);

  template <class DerivedIn, class DerivedOut, int BlockRows, int BlockCols, bool InnerPanel> void so2Exp(Eigen::MatrixBase<DerivedIn> const& delta, Eigen::Block<DerivedOut, BlockRows, BlockCols, InnerPanel> res){

    using S = typename DerivedIn::Scalar;

    S const cd{ std::cos(delta(0)) };
    S const sd{ std::sin(delta(0)) };

    res << cd, -sd,
      sd, cd;    
  }

  template <class S, class DerivedOut, class RetType = void> 
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta, Eigen::MatrixBase<DerivedOut>& res);

  template <class S, class DerivedOut, int BlockRows, int BlockCols, bool InnerPanel, class RetType = void> 
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta, Eigen::Block<DerivedOut, BlockRows, BlockCols, InnerPanel> res){
    S const cd{ std::cos(delta) };
    S const sd{ std::sin(delta) };

    res << cd, -sd,
      sd, cd;    
  }

  // @brief compute the exp of a so2 algebra to a SO2
  // @tparam S scalar type
  // @param delta a scalar vector representing the so2 algebra vector
  // @return the rotation matrix in SO2 from the algebra vector in so2

  template <class S, class RetType = SO2<S>>
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta);


}