#pragma once

#include <Eigen/Core>

#include <csCore/BasicTypes.hpp>

namespace csLie {

  // S2, a 3D unit vector
  // note, there are no guarantees that S2 is unitary, is up to the user to keep it coherent
  template <class S> using S2 = Eigen::Matrix<S, 3, 1>;
  using S2d = S2<float64_t>;
  using S2f = S2<float32_t>;
  
  // S2 algebra, a 2D vector
  template <class S> using S2Alg = Eigen::Matrix<S, 2, 1>;
  using S2Algd = S2Alg<float64_t>;
  using S2Algf = S2Alg<float32_t>;

  // S2 jacobian for oplus operation, a 3x2 matrix
  // vp = oplusS2(S2 v, S2Alg delta)
  template <class S> using S2OPlusJac = Eigen::Matrix<S, 3, 2>;
  using S2OPlusJacd = S2OPlusJac<float64_t>;
  using S2OPlusJacf = S2OPlusJac<float32_t>;


  // SO3 = 3D rotation matrix, a 3x3 matrix orthogonal with det = 1
  // note, there are no guarantees that SO3 is a proper rotation matrix, is up to the user to keep it coherent
  template <class S> using SO3 = Eigen::Matrix<S, 3, 3>;
  using SO3d = SO3<float64_t>;
  using SO3f = SO3<float32_t>;

}