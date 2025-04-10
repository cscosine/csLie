#pragma once

#include <Eigen/Core>
#include <csLie/Tagged.hpp>

namespace csLie {

template <class S>
using SE3Alg = Tagged<S, 6, 1, struct TagSE3Alg>;

using SE3Algd = SE3Alg<double>;
using SE3Algf = SE3Alg<float>;

}  // namespace csLie
