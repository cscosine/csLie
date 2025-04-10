#pragma once
#include <Eigen/Core>
#include <csLie/SE3.hpp>
#include <csLie/Tagged.hpp>

namespace csLie {

struct TagSE3Params {};

// TODO, do not extend, compose instead
template <class S>
class SE3Params : public Tagged<S, 6, 1, TagSE3Params> {
 public:
  SE3Params() {}

  template <class Derived>
  explicit SE3Params(const Eigen::EigenBase<Derived>& d)
      : Tagged<S, 6, 1, TagSE3Params>(d) {}

  explicit SE3Params(S x, S y, S z, S yaw, S pitch, S roll);

  Eigen::Block<Eigen::Matrix<S, 6, 1>, 3, 1> translation() {
    return this->template head<3>();
  }

  const Eigen::Block<const Eigen::Matrix<S, 6, 1>, 3, 1> translation() const {
    return this->template head<3>();
  }

  Eigen::Block<Eigen::Matrix<S, 6, 1>, 3, 1> eulerAngles() {
    return this->template tail<3>();
  }

  const Eigen::Block<const Eigen::Matrix<S, 6, 1>, 3, 1> eulerAngles() const {
    return this->template tail<3>();
  }
};

template <typename S>
Pose3D<S> params2Pose(const SE3Params<S>& g);

template <typename S>
SE3Params<S> pose2Params(const Pose3D<S>& g);

using SE3Paramsd = SE3Params<double>;
using SE3Paramsf = SE3Params<float>;

}  // namespace csLie
