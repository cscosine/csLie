#pragma once

#include <Eigen/Core>

#include <csCore/BasicTypes.hpp>

namespace csLie {

  template <class S>
  using S2 = Eigen::Matrix<S, 3, 1>;

  using S2d = S2<float64_t>;
  using S2f = S2<float32_t>;

  template <class S>
  Eigen::Matrix<S, 3, 3> s2toRotationMatrix(const S2<S> & vec);

}