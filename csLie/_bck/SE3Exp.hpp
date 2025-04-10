#pragma once
#include <csLie/SE3.hpp>
#include <csLie/SE3Alg.hpp>

namespace csLie {

template <typename S>
Pose3D<S> exp(const SE3Alg<S>& g);

}
