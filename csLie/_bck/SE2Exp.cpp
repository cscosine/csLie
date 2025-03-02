#include <csLie/SE2Exp.hpp>


namespace csLie {

  template<typename S>
  Pose2D<S> exp(const SE2Alg<S>& g) {


    S cosTh = std::cos(g(2));
    S sinTh = std::sin(g(2));
    S sinTh_Th, m1cosTh_Th;
    if (std::abs(g(2)) < Eigen::NumTraits<S>::epsilon()) {
      sinTh_Th = 1;
      m1cosTh_Th = 0;
    }
    else {
      sinTh_Th = sinTh / g(2);
      m1cosTh_Th = (1 - cosTh) / g(2);
    }

    Pose2D<S> ret;
    ret.linear()(0, 0) = cosTh;
    ret.linear()(0, 1) = -sinTh;
    ret.linear()(1, 0) = sinTh;
    ret.linear()(1, 1) = cosTh;

    ret.translation().x() = g.x() * sinTh_Th - g.y() * m1cosTh_Th;
    ret.translation().y() = g.x() * m1cosTh_Th + g.y() * sinTh_Th;
    return ret;


  }

  template Pose2D<double> exp<double>(const SE2Alg<double>& g);
  template Pose2D<float> exp<float>(const SE2Alg<float>& g);

}

