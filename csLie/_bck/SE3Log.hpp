#pragma once
#include <csLie/SE3.hpp>
#include <csLie/SE3Alg.hpp>

namespace csLie {

  template<typename S>
  SE3Alg<S> log(const Pose3D<S>& g);

}
