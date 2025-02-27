#pragma once
#include "SE3.h"
#include "SE3Alg.h"

namespace csLie {

  template<typename S>
  SE3Alg<S> log(const Pose3D<S>& g);

}
