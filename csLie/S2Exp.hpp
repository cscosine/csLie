#pragma once

#include <Eigen/Core>

#include <csLie/S2.hpp>
#include <csLie/S2Alg.hpp>

namespace csLie {

  template<class S>
  S2<S> exp(S2Alg<S> in);

  template<class S>
  Eigen::Matrix<S, 3, 2> d_Rv_expd_d_expd0(const S2<S>& v);

}