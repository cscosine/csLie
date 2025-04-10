#pragma once

#include <Eigen/Core>
#include <csCore/core/BasicTypes.hpp>

namespace csLie {

/// @brief a set of constants used in the library
/// @tparam S the scalar type
template <typename S>
class Constants {
 public:
  /// @brief the machine epsilon, i.e. the minium number that makes (1 + epsilon
  /// == 1) true
  static constexpr S epsilon{Eigen::NumTraits<S>::epsilon()};

  /// @brief g a weak epsilon value. It is mainly used as a default value for
  /// comparison operators.
  static constexpr S dummy_precision{Eigen::NumTraits<S>::dummy_precision()};

  /// @brief pi, what else
  static constexpr S pi{static_cast<S>(EIGEN_PI)};
};

}  // namespace csLie
