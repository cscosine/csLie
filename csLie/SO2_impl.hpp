#pragma once

#include <csLie/SO2.hpp>
#include <csLie/SO2.hpp>
#include <csCore/eigenUtils/MatrixShape.hpp>
#include <csLie/csLieAssert.hpp>

#include <Eigen/Dense>

namespace csLie {

  template <class Derived> 
  SO2<typename Derived::Scalar> so2Exp(Eigen::MatrixBase<Derived> const& delta) {
    using S = typename Derived::Scalar;
    using M11 = Eigen::Matrix<S, 1, 1>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(delta, M11()), "delta required to be 1x1");
    return so2Exp(delta(0,0));
  }

  template <class S, class RetType>
  std::enable_if_t<std::is_floating_point<S>::value, RetType>
  so2Exp(S delta) {
    SO2<S> res{};
    S const cd{ std::cos(delta) };
    S const sd{ std::sin(delta) };

    res(0,0) = cd;
    res(0,1) = -sd;
    res(1,0) = sd;
    res(1,1) = cd;
    return res;
  }

  template <class Derived> 
  bool isValidSO2RotationMatrix(Eigen::MatrixBase<Derived> const& rot, typename Derived::Scalar precision)
  {
    using S = typename Derived::Scalar;

    using M22 = Eigen::Matrix<S, 2, 2>;
    CSLIE_CONTRACT_DEBUG(precision >= static_cast<S>(0.0), "precision needs to be non negative");
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(rot, M22()), "rot required to be 2x2");

    S const det{rot.determinant()};
    if(std::abs(det - static_cast<S>(1.0)) >= precision) {
      return false;
    }
    
    if(!(rot * rot.transpose() - Eigen::Matrix2<S>::Identity()).isZero(precision)) {
      return false;
    }

    if(!(rot.transpose() * rot - Eigen::Matrix2<S>::Identity()).isZero(precision)) {
      return false;
    }

    return true;
  }

  template <class Derived>
  typename Derived::Scalar so2Log(Eigen::MatrixBase<Derived> const& rot) 
  {
    using S = typename Derived::Scalar;

    using M22 = Eigen::Matrix<S, 2, 2>;
    CSLIE_CONTRACT_DEBUG(csCore::eigenUtils::isSameSizeStaticOrDynamic(rot, M22()), "rot required to be 2x2");
    CSLIE_CONTRACT_PARANOID(isValidSO2RotationMatrix(rot), "rot is not a valid rotation matrix");
    
    // atan2 (sin,cos), average of values 
    return std::atan2(rot(1,0), rot(0,0));
  }


}