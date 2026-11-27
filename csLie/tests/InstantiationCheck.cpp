// Compiles every public function of csLie in float32_t and float64_t under strict warnings (cscmake_strict_warnings,
// set on this target in CMakeLists.txt). csLie is header-only: a template body is compiled, and warned about, only
// where it is instantiated, and the tests build at a lower warning level. Nothing here runs.
// A new public function is added here.
#include <csLie/S2_impl.hpp>
#include <csLie/SE2_impl.hpp>
#include <csLie/SE3_impl.hpp>
#include <csLie/SO2_impl.hpp>
#include <csLie/SO3_impl.hpp>
#include <csLie/Types.hpp>
#include <csLie/Utils_impl.hpp>

namespace csLie {

template <class S> void instantiateAll() {
  SO2Alg<S> const so2Alg = SO2Alg<S>::Zero();
  SO2<S> const so2 = SO2<S>::Identity();
  SE2Alg<S> const se2Alg = SE2Alg<S>::Zero();
  SE2<S> const se2 = SE2<S>::Identity();
  SO3Alg<S> const so3Alg = SO3Alg<S>::Zero();
  SO3<S> const so3 = SO3<S>::Identity();
  SE3Alg<S> const se3Alg = SE3Alg<S>::Zero();
  SE3<S> const se3 = SE3<S>::Identity();
  S2<S> const s2 = S2<S>::UnitX();
  S2Alg<S> const s2Alg = S2Alg<S>::Zero();

  // SO2
  static_cast<void>(so2Exp(so2Alg));
  static_cast<void>(so2Exp(S(0)));
  static_cast<void>(so2Log(so2));
  // SE2
  static_cast<void>(se2Exp(se2Alg));
  static_cast<void>(se2Log(se2.matrix()));
  static_cast<void>(se2Log(se2));
  static_cast<void>(xyYToSe2(se2Alg));
  static_cast<void>(se2ToxyY(se2.matrix()));
  static_cast<void>(se2ToxyY(se2));
  // SO3
  static_cast<void>(so3Exp(so3Alg));
  static_cast<void>(so3Log(so3));
  static_cast<void>(so3J(so3Alg));
  static_cast<void>(so3Jinv(so3Alg));
  static_cast<void>(so3ToRPY(so3));
  static_cast<void>(rpyToSo3(so3Alg));
  static_cast<void>(skewSymm(so3Alg));
  // SE3
  static_cast<void>(se3Exp(se3Alg));
  static_cast<void>(se3Log(se3.matrix()));
  static_cast<void>(se3Log(se3));
  static_cast<void>(se3Adjoint(se3.matrix()));
  static_cast<void>(se3Adjoint(se3));
  static_cast<void>(se3JacobianLeftZeroPerturbationPointTransformation(so3Alg));
  // S2
  static_cast<void>(s2Exp(s2Alg));
  static_cast<void>(s2toZAlignedSO3(s2));
  static_cast<void>(s2oplus(s2, s2Alg));
  static_cast<void>(s2_d_oplus_d_delta_0(s2));
  // Utils
  static_cast<void>(normalize3dVectorWithJacobians(s2));
}

template void instantiateAll<float32_t>();
template void instantiateAll<float64_t>();

} // namespace csLie
