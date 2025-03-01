#include <csLie/S2.hpp>

namespace csLie {

  template <class S>
  Eigen::Matrix<S, 3, 3> s2toRotationMatrix(const S2<S>& vec) {

    S const alpha = std::atan2(vec.z(), vec.y());
    S const r = std::sqrt(vec.z() * vec.z() + vec.y() * vec.y());
    S const c = cos(alpha);
    S const s = sin(alpha);

    Eigen::Matrix<S, 3, 3> R;
    R.col(0) = vec;
    R.col(1) = Eigen::Matrix<S, 3, 1>(-r, vec.x() * c, vec.x() * s);
    R.col(2) = Eigen::Matrix<S, 3, 1>(0, -s, c);

    return R;
  }

  template Eigen::Matrix<float64_t, 3, 3> s2toRotationMatrix(const S2<float64_t>& vec);
  template Eigen::Matrix<float32_t, 3, 3> s2toRotationMatrix(const S2<float32_t>& vec);


}