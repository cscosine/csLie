#include <csLie/SE3Params.hpp>
#include <csLie/SO3Params.hpp>

namespace csLie {

template <typename S>
SE3Params<S>::SE3Params(S x, S y, S z, S yaw, S pitch, S roll) {
  (*this)(0) = x;
  (*this)(1) = y;
  (*this)(2) = z;
  (*this)(3) = yaw;
  (*this)(4) = pitch;
  (*this)(5) = roll;
}

//-------------------------------------------

template <typename S>
Pose3D<S> params2Pose(const SE3Params<S>& g) {
  Pose3D<S> T = Pose3D<S>::Identity();
  T.translation() = g.translation();
  T.linear() = params2Pose(SO3Params<S>{g.eulerAngles()});
  return T;
}

template <typename S>
SE3Params<S> pose2Params(const Pose3D<S>& g) {
  SE3Params<S> ret;
  ret.translation() = g.translation();
  ret.eulerAngles() = pose2Params(Rot3D<S>{g.linear()});
  return ret;
}

template class SE3Params<double>;
template class SE3Params<float>;

template Pose3D<double> params2Pose<double>(const SE3Params<double>& g);
template Pose3D<float> params2Pose<float>(const SE3Params<float>& g);

template SE3Params<double> pose2Params<double>(const Pose3D<double>& g);
template SE3Params<float> pose2Params<float>(const Pose3D<float>& g);

}  // namespace csLie
