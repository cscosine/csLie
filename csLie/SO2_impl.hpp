#pragma once

#include <csLie/SO2.hpp>
#include <csLie/SO2.hpp>
#include <csCore/eigenUtils/MatrixShape.hpp>
#include <csLie/csLieAssert.hpp>

namespace csLie {

  namespace _private {
    template<class S, class O>
    void so2ExpImpl(S delta, O res) {
      S const cd{ std::cos(delta) };
      S const sd{ std::sin(delta) };
  
      res(0,0) = cd;
      res(0,1) = -sd;
      res(1,0) = sd;
      res(1,1) = cd;
  
    }
  }

  template <class Derived> 
  SO2<typename Derived::Scalar> so2Exp(Eigen::MatrixBase<Derived> const& delta) {
    using S = typename Derived::Scalar;
    using M11 = Eigen::Matrix<S, 1, 1>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, M11()), "delta required to be 1x1");
    
    SO2<S> r{};
    so2Exp(delta, r);
    return r;
  }

  template <class DerivedIn, class DerivedOut> 
  void so2Exp(Eigen::MatrixBase<DerivedIn> const& delta, Eigen::MatrixBase<DerivedOut>& res)
  {
    using S = typename DerivedIn::Scalar;
    using M11 = Eigen::Matrix<S, 1, 1>;
    using M22 = Eigen::Matrix<S, 2, 2>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, M11()), "delta required to be 1x1");
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(res, M22()), "output required to be 2x2");


    _private::so2ExpImpl<S,Eigen::MatrixBase<DerivedOut>&>(delta(0,0), res);
  }

  template <class DerivedIn, class DerivedOut, int BlockRows, int BlockCols, bool InnerPanel> 
  void so2Exp(Eigen::MatrixBase<DerivedIn> const& delta, Eigen::Block<DerivedOut, BlockRows, BlockCols, InnerPanel> res){

    using S = typename DerivedIn::Scalar;

    using M11 = Eigen::Matrix<S, 1, 1>;
    using M22 = Eigen::Matrix<S, 2, 2>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, M11()), "delta required to be 1x1");
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(res, M22()), "output required to be 2x2");

    _private::so2ExpImpl<S,Eigen::Block<DerivedOut, BlockRows, BlockCols, InnerPanel> >(delta(0,0), res);

  }

  template <class S, class RetType>
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta) {
    return so2Exp(SO2Alg<S>(delta));
  }

  template <class S, class DerivedOut, class RetType>
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
    so2Exp(S delta, Eigen::MatrixBase<DerivedOut>& res)
  {
    so2Exp(SO2Alg<S>(delta), res);
  }

  template <class S, class DerivedOut, int BlockRows, int BlockCols, bool InnerPanel, class RetType> 
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta, Eigen::Block<DerivedOut, BlockRows, BlockCols, InnerPanel> res){
    using M22 = Eigen::Matrix<S, 2, 2>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(res, M22()), "output required to be 2x2");

    _private::so2ExpImpl<S,Eigen::Block<DerivedOut, BlockRows, BlockCols, InnerPanel> >(delta, res);
  }

  template <class DerivedIn, int BlockRows, int BlockCols, bool InnerPanel> 
  SO2<typename DerivedIn::Scalar> so2Exp(Eigen::Block<DerivedIn, BlockRows, BlockCols, InnerPanel> const& delta) {
    using S = typename DerivedIn::Scalar;

    using M11 = Eigen::Matrix<S, 1, 1>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, M11()), "delta required to be 1x1");

    SO2<S> r{};
    so2Exp(SO2Alg<S>(delta(0,0)), r);
    return r;
  }

  template <class DerivedIn, int BlockRows, int BlockCols, bool InnerPanel, class DerivedOut> 
  void so2Exp(Eigen::Block<DerivedIn, BlockRows, BlockCols, InnerPanel> const& delta, Eigen::MatrixBase<DerivedOut> & res)
  {
    using Sin = typename DerivedIn::Scalar;
    using Sout = typename DerivedIn::Scalar;
    static_assert(std::is_same<Sin, Sout>::value, "so2Exp, input and output types should be the same");

    using S = Sin;

    using M11 = Eigen::Matrix<S, 1, 1>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, M11()), "delta required to be 1x1");
    
    so2Exp(SO2Alg<S>(delta(0,0)), res);

  }

  template <class DerivedIn, int BlockRowsIn, int BlockColsIn, bool InnerPanelIn, class DerivedOut, int BlockRowsOut, int BlockColsOut, bool InnerPanelOut> 
  void so2Exp(Eigen::Block<DerivedIn, BlockRowsIn, BlockColsIn, InnerPanelIn> const& delta, Eigen::Block<DerivedOut, BlockRowsOut, BlockColsOut, InnerPanelOut> res)
  {
    using Sin = typename DerivedIn::Scalar;
    using Sout = typename DerivedIn::Scalar;
    static_assert(std::is_same<Sin, Sout>::value, "so2Exp, input and output types should be the same");

    using S = Sin;

    using M11 = Eigen::Matrix<S, 1, 1>;
    using M22 = Eigen::Matrix<S, 2, 2>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, M11()), "delta required to be 1x1");
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(res, M22()), "output required to be 2x2");
    
    so2Exp(SO2Alg<S>(delta(0,0)), res);

  }

}