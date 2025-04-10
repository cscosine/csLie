#pragma once
#include <Eigen/Core>
#include <csLie/Tagged.hpp>

namespace csLie {

template <class S>
using SO3Alg = Tagged<S, 3, 1, struct TagSO3Alg>;

using SO3Algd = SO3Alg<double>;
using SO3Algf = SO3Alg<float>;

}  // namespace csLie
