#pragma once
#include "SE3.h"
#include "SE3Alg.h"

namespace csLie {

  template<typename S>
  Pose3D<S> exp(const SE3Alg<S>& g);

}
