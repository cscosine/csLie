#pragma once

#include <Eigen/Core>
#include <csCore/core/BasicTypes.hpp>

namespace csLie {

// ------------ S2 ------------

/// @brief S2, a 3D unit vector
/// note, there are no guarantees that S2 is unitary, is up to the user to keep
/// it coherent
/// @tparam S the scalar type
template <class S>
using S2 = Eigen::Matrix<S, 3, 1>;

/// @brief convenient typedef for float64_t
using S2d = S2<float64_t>;
/// @brief convenient typedef for float32_t
using S2f = S2<float32_t>;

/// @brief S2 algebra, a 2D vector
/// @tparam S the scalar type
template <class S>
using S2Alg = Eigen::Matrix<S, 2, 1>;
/// @brief convenient typedef for float64_t
using S2Algd = S2Alg<float64_t>;
/// @brief convenient typedef for float32_t
using S2Algf = S2Alg<float32_t>;

/// @brief the matrix storing the Jacobian of s2 oplus wrt delta evaluated in
/// delta = 0
/// @tparam S the scalar type
template <class S>
using s2oplusJac = Eigen::Matrix<S, 3, 2>;
/// @brief convenient typedef for float64_t
using s2oplusJacd = s2oplusJac<float64_t>;
/// @brief convenient typedef for float32_t
using s2oplusJacf = s2oplusJac<float32_t>;

// ------------ SO2 ------------

/// @brief SO2 = 2D rotation matrix, a 2x2 matrix orthogonal with det = 1
/// note, there are no guarantees that SO2 is a proper rotation matrix, is up to
/// the user to keep it coherent
/// @tparam S the scalar type
template <class S>
using SO2 = Eigen::Matrix<S, 2, 2>;
/// @brief convenient typedef for float64_t
using SO2d = SO2<float64_t>;
/// @brief convenient typedef for float32_t
using SO2f = SO2<float32_t>;

/// @brief S2 algebra, a 1D vector
/// note: despite this can be represented with a basic type, we prefer for
/// coherence to use a 1D vector
/// @tparam S the scalar type
template <class S>
using SO2Alg = Eigen::Matrix<S, 1, 1>;
/// @brief convenient typedef for float64_t
using SO2Algd = SO2Alg<float64_t>;
/// @brief convenient typedef for float32_t
using SO2Algf = SO2Alg<float32_t>;

// ------------ SE2 ------------

template <class S>
using SE2 = Eigen::Transform<S, 2, Eigen::Isometry>;
/// @brief convenient typedef for float64_t
using SE2d = SE2<float64_t>;
/// @brief convenient typedef for float32_t
using SE2f = SE2<float32_t>;

/// @brief SE2 algebra, a 3D vector
/// @tparam S the scalar type
template <class S>
using SE2Alg = Eigen::Matrix<S, 3, 1>;
/// @brief convenient typedef for float64_t
using SE2Algd = SE2Alg<float64_t>;
/// @brief convenient typedef for float32_t
using SE2Algf = SE2Alg<float32_t>;

// ------------ SO3 ------------

/// @brief SO3 = 3D rotation matrix, a 3x3 matrix orthogonal with det = 1
/// note, there are no guarantees that SO3 is a proper rotation matrix, is up to
/// the user to keep it coherent
/// @tparam S the scalar type
template <class S>
using SO3 = Eigen::Matrix<S, 3, 3>;
/// @brief convenient typedef for float64_t
using SO3d = SO3<float64_t>;
/// @brief convenient typedef for float32_t
using SO3f = SO3<float32_t>;

}  // namespace csLie
