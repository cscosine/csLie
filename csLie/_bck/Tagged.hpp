#pragma once
#include <Eigen/Core>

namespace csLie {

template <class T, int DimRows, int DimCols, class Tag>
class Tagged : public Eigen::Matrix<T, DimRows, DimCols> {
  //
};

}  // namespace csLie
