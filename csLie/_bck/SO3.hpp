#pragma once
#include <Eigen/Core>
#include <csLie/Tagged.hpp>

namespace csLie {

template <class S>
using Rot3D = Tagged<S, 3, 3, struct TagSO3>;

using Rot3Dd = Rot3D<double>;
using Rot3Df = Rot3D<float>;

}  // namespace csLie
