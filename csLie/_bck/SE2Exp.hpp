#pragma once
#include <csLie/SE2.hpp>
#include <csLie/SE2Alg.hpp>

namespace csLie {

template <typename S>
Pose2D<S> exp(const SE2Alg<S>& g);

}
