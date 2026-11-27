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
template <class S> using S2 = Eigen::Matrix<S, 3, 1>;

/// @brief convenient typedef for float64_t
using S2d = S2<float64_t>;
/// @brief convenient typedef for float32_t
using S2f = S2<float32_t>;

/// @brief S2 tangent vector, 2D: coordinates in the tangent basis of s2Frame (S2 is not a group, no Lie algebra)
/// @tparam S the scalar type
template <class S> using S2Tangent = Eigen::Matrix<S, 2, 1>;
/// @brief convenient typedef for float64_t
using S2Tangentd = S2Tangent<float64_t>;
/// @brief convenient typedef for float32_t
using S2Tangentf = S2Tangent<float32_t>;

// ------------ SO2 ------------

/// @brief SO2 = 2D rotation matrix, a 2x2 matrix orthogonal with det = 1
/// note, there are no guarantees that SO2 is a proper rotation matrix, is up to
/// the user to keep it coherent
/// @tparam S the scalar type
template <class S> using SO2 = Eigen::Matrix<S, 2, 2>;
/// @brief convenient typedef for float64_t
using SO2d = SO2<float64_t>;
/// @brief convenient typedef for float32_t
using SO2f = SO2<float32_t>;

/// @brief SO2 tangent vector (Lie algebra coordinates), 1D
/// note: despite this can be represented with a basic type, we prefer for
/// coherence to use a 1D vector
/// @tparam S the scalar type
template <class S> using SO2Tangent = Eigen::Matrix<S, 1, 1>;
/// @brief convenient typedef for float64_t
using SO2Tangentd = SO2Tangent<float64_t>;
/// @brief convenient typedef for float32_t
using SO2Tangentf = SO2Tangent<float32_t>;

// ------------ SE2 ------------

/// @brief SE2 = 2D transformation matrix, inside it is a 3x3 matrix, top left
/// 2x2 rotation, last column translation, last row 0 0 1 note, there are no
/// guarantees that SE2 is a proper rotation matrix, is up to the user to keep
/// it coherent
/// @tparam S the scalar type
template <class S> using SE2 = Eigen::Transform<S, 2, Eigen::Isometry>;
/// @brief convenient typedef for float64_t
using SE2d = SE2<float64_t>;
/// @brief convenient typedef for float32_t
using SE2f = SE2<float32_t>;

/// @brief SE2 in matrix form
template <class S> using SE2Matrix = Eigen::Matrix<S, 3, 3>;
/// @brief convenient typedef for float64_t
using SE2Matrixd = SE2Matrix<float64_t>;
/// @brief convenient typedef for float32_t
using SE2Matrixf = SE2Matrix<float32_t>;

/// @brief SE2 tangent vector (Lie algebra coordinates), 3D
/// @tparam S the scalar type
template <class S> using SE2Tangent = Eigen::Matrix<S, 3, 1>;
/// @brief convenient typedef for float64_t
using SE2Tangentd = SE2Tangent<float64_t>;
/// @brief convenient typedef for float32_t
using SE2Tangentf = SE2Tangent<float32_t>;

// ------------ SO3 ------------

/// @brief SO3 = 3D rotation matrix, a 3x3 matrix orthogonal with det = 1
/// note, there are no guarantees that SO3 is a proper rotation matrix, is up to
/// the user to keep it coherent
/// @tparam S the scalar type
template <class S> using SO3 = Eigen::Matrix<S, 3, 3>;
/// @brief convenient typedef for float64_t
using SO3d = SO3<float64_t>;
/// @brief convenient typedef for float32_t
using SO3f = SO3<float32_t>;

/// @brief SO3 tangent vector (Lie algebra coordinates), 3D
/// @tparam S the scalar type
template <class S> using SO3Tangent = Eigen::Matrix<S, 3, 1>;
/// @brief convenient typedef for float64_t
using SO3Tangentd = SO3Tangent<float64_t>;
/// @brief convenient typedef for float32_t
using SO3Tangentf = SO3Tangent<float32_t>;

// ------------ SE3 ------------

/// @brief SE3 = 3D transformation matrix, a 4x4 matrix with top left rotation
/// and translation vector in last column note, there are no guarantees that SE3
/// is a proper rotation matrix, is up to the user to keep it coherent
/// @tparam S the scalar type
template <class S> using SE3 = Eigen::Transform<S, 3, Eigen::Isometry>;
/// @brief convenient typedef for float64_t
using SE3d = SE3<float64_t>;
/// @brief convenient typedef for float32_t
using SE3f = SE3<float32_t>;

/// @brief SE3 in matrix form
template <class S> using SE3Matrix = Eigen::Matrix<S, 4, 4>;
/// @brief convenient typedef for float64_t
using SE3Matrixd = SE3Matrix<float64_t>;
/// @brief convenient typedef for float32_t
using SE3Matrixf = SE3Matrix<float32_t>;

/// @brief SE3 tangent vector (Lie algebra coordinates), 6D
/// @tparam S the scalar type
template <class S> using SE3Tangent = Eigen::Matrix<S, 6, 1>;
/// @brief convenient typedef for float64_t
using SE3Tangentd = SE3Tangent<float64_t>;
/// @brief convenient typedef for float32_t
using SE3Tangentf = SE3Tangent<float32_t>;

/// @brief SE3 adjoint matrix, a 6x6 matrix
/// @tparam S the scalar type
template <class S> using SE3Adjoint = Eigen::Matrix<S, 6, 6>;
/// @brief convenient typedef for float64_t
using SE3Adjointd = SE3Adjoint<float64_t>;
/// @brief convenient typedef for float32_t
using SE3Adjointf = SE3Adjoint<float32_t>;

// ------------ Points ------------
/// @brief 3d point, a 3x1 vector
/// @tparam S the scalar type
template <class S> using Point3d = Eigen::Matrix<S, 3, 1>;

/// @brief 2d point, a 2x1 vector
/// @tparam S the scalar type
template <class S> using Point2d = Eigen::Matrix<S, 2, 1>;

// ------------ Utils ------------

/// @brief a struct containing a unit vector and the jacobian computed to
/// normalize the input vector, used as return type of
/// normalizeWithJacobian
template <class S> struct NormalizedWithJacobian {
  Eigen::Matrix<S, 3, 1> unit;
  Eigen::Matrix<S, 3, 3> jacobian;
};
/// @brief convenient typedef for float64_t
using NormalizedWithJacobiand = NormalizedWithJacobian<float64_t>;
/// @brief convenient typedef for float32_t
using NormalizedWithJacobianf = NormalizedWithJacobian<float32_t>;

} // namespace csLie
