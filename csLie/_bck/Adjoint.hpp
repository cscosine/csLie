#pragma once
#include <Eigen/Core>
#include <csLie/SE3.hpp>


namespace csLie {

  template<class S>
  Eigen::Matrix<S, 6, 6> adjoint(const Pose3D<S>& p);

}

