#pragma once

#include <Eigen/Core>

namespace csLie {
  // perform a set of checks on a rotation matrix to verify its properties:
  // - det = 1
  // - cols and rows pairwise orthogonal
  // - col and row norm = 1
  // - R^T * R = I
  // - cols are right handed 
  template <class S> void checkValidRotationMatrix(Eigen::Matrix3<S> const& rot);
}