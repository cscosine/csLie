#pragma once
#include <csLie/SO3.hpp>
#include <csLie/SO3Alg.hpp>

namespace csLie {

  template<typename S>
  SO3Alg<S> log(const Rot3D<S>& g);

}
