#pragma once
#include <csLie/SO3.hpp>
#include <csLie/SO3Alg.hpp>

namespace csLie {

  template<typename S>
  Rot3D<S> exp(const SO3Alg<S>& g);

}
