#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
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

/// @brief SE2 = 2D transformation matrix, inside it is a 3x3 matrix, top left
/// 2x2 rotation, last column translation, last row 0 0 1 note, there are no
/// guarantees that SE2 is a proper rotation matrix, is up to the user to keep
/// it coherent
/// @tparam S the scalar type
template <class S>
using SE2 = Eigen::Transform<S, 2, Eigen::Isometry>;
/// @brief convenient typedef for float64_t
using SE2d = SE2<float64_t>;
/// @brief convenient typedef for float32_t
using SE2f = SE2<float32_t>;

/// @brief SE2 in matrix form
template <class S>
using SE2Matrix = Eigen::Matrix<S, 3, 3>;
/// @brief convenient typedef for float64_t
using SE2Matrixd = SE2Matrix<float64_t>;
/// @brief convenient typedef for float32_t
using SE2Matrixf = SE2Matrix<float32_t>;

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

/// @brief SO3 algebra, a 3D vector
/// @tparam S the scalar type
template <class S>
using SO3Alg = Eigen::Matrix<S, 3, 1>;
/// @brief convenient typedef for float64_t
using SO3Algd = SO3Alg<float64_t>;
/// @brief convenient typedef for float32_t
using SO3Algf = SO3Alg<float32_t>;

// ------------ SE3 ------------

/// @brief SE3 = 3D transformation matrix, a 4x4 matrix with top left rotation
/// and translation vector in last column note, there are no guarantees that SE3
/// is a proper rotation matrix, is up to the user to keep it coherent
/// @tparam S the scalar type
template <class S>
using SE3 = Eigen::Transform<S, 3, Eigen::Isometry>;
/// @brief convenient typedef for float64_t
using SE3d = SE3<float64_t>;
/// @brief convenient typedef for float32_t
using SE3f = SE3<float32_t>;

/// @brief SE3 in matrix form
template <class S>
using SE3Matrix = Eigen::Matrix<S, 4, 4>;
/// @brief convenient typedef for float64_t
using SE3Matrixd = SE3Matrix<float64_t>;
/// @brief convenient typedef for float32_t
using SE3Matrixf = SE3Matrix<float32_t>;

/// @brief SE3 algebra, a 6D vector
/// @tparam S the scalar type
template <class S>
using SE3Alg = Eigen::Matrix<S, 6, 1>;
/// @brief convenient typedef for float64_t
using SE3Algd = SE3Alg<float64_t>;
/// @brief convenient typedef for float32_t
using SE3Algf = SE3Alg<float32_t>;

/// @brief SE3 andjoint matrix, a 6x6 vector
/// @tparam S the scalar type
template <class S>
using SE3Adj = Eigen::Matrix<S, 6, 6>;
/// @brief convenient typedef for float64_t
using SE3Adjd = SE3Adj<float64_t>;
/// @brief convenient typedef for float32_t
using SE3Adjf = SE3Adj<float32_t>;

// ------------ Utils ------------

/// @brief a struct containing a unit vector and the jacobian computed to
/// normalize the input vector, used as return type of
/// normalize3dVectorWithJacobians
template <class S>
struct UnitVector3dWithJacobians {
  Eigen::Matrix<S, 3, 1> vo;
  Eigen::Matrix<S, 3, 3> d_vo_d_v;
};
/// @brief convenient typedef for float64_t
using UnitVector3dWithJacobiansd = UnitVector3dWithJacobians<float64_t>;
/// @brief convenient typedef for float32_t
using UnitVector3dWithJacobiansf = UnitVector3dWithJacobians<float32_t>;

}  // namespace csLie
