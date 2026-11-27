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
  SO2Tangent<S> const so2Tangent = SO2Tangent<S>::Zero();
  SO2<S> const so2 = SO2<S>::Identity();
  SE2Tangent<S> const se2Tangent = SE2Tangent<S>::Zero();
  SE2<S> const se2 = SE2<S>::Identity();
  SO3Tangent<S> const so3Tangent = SO3Tangent<S>::Zero();
  SO3<S> const so3 = SO3<S>::Identity();
  SE3Tangent<S> const se3Tangent = SE3Tangent<S>::Zero();
  SE3<S> const se3 = SE3<S>::Identity();
  S2<S> const s2 = S2<S>::UnitX();
  S2Tangent<S> const s2Tangent = S2Tangent<S>::Zero();

  // SO2
  static_cast<void>(so2Exp(so2Tangent));
  static_cast<void>(so2Exp(S(0)));
  static_cast<void>(so2Log(so2));
  // SE2
  static_cast<void>(se2Exp(se2Tangent));
  static_cast<void>(se2Log(se2.matrix()));
  static_cast<void>(se2Log(se2));
  static_cast<void>(xyYawToSe2(se2Tangent));
  static_cast<void>(xyYawToSe2(Eigen::Matrix<S, 2, 1>::Zero().eval(), S(0)));
  static_cast<void>(xyYawToSe2(S(0), S(0), S(0)));
  static_cast<void>(se2ToXyYaw(se2.matrix()));
  static_cast<void>(se2ToXyYaw(se2));
  // SO3
  static_cast<void>(so3Exp(so3Tangent));
  static_cast<void>(so3Log(so3));
  static_cast<void>(so3LeftJacobian(so3Tangent));
  static_cast<void>(so3LeftJacobianInverse(so3Tangent));
  static_cast<void>(so3RightJacobian(so3Tangent));
  static_cast<void>(so3RightJacobianInverse(so3Tangent));
  static_cast<void>(so3OplusLeft(so3, so3Tangent));
  static_cast<void>(so3OplusRight(so3, so3Tangent));
  static_cast<void>(vee(so3));
  static_cast<void>(so3ToRollPitchYaw(so3));
  static_cast<void>(so3ToRollPitchYaw<GimbalPolicy::kYawZero>(so3));
  static_cast<void>(so3ToRollPitchYawChecked(so3));
  static_cast<void>(rollPitchYawCheckedGimbalTolerance<S>());
  static_cast<void>(rollPitchYawToSo3(so3Tangent));
  static_cast<void>(hat(so3Tangent));
  // SE3
  static_cast<void>(se3Exp(se3Tangent));
  static_cast<void>(se3Log(se3.matrix()));
  static_cast<void>(se3Log(se3));
  static_cast<void>(se3Adjoint(se3.matrix()));
  static_cast<void>(se3Adjoint(se3));
  static_cast<void>(se3PointJacobianLeft(so3Tangent));
  static_cast<void>(se3PointJacobianRight(so3, so3Tangent));
  static_cast<void>(se3OplusLeft(se3, se3Tangent));
  static_cast<void>(se3OplusRight(se3, se3Tangent));
  // S2
  static_cast<void>(s2Exp(s2Tangent));
  static_cast<void>(s2Frame(s2));
  static_cast<void>(s2Oplus(s2, s2Tangent));
  static_cast<void>(s2OplusJacobianAtZero(s2));
  // Utils
  static_cast<void>(normalizeWithJacobian(s2));
}

template void instantiateAll<float32_t>();
template void instantiateAll<float64_t>();

} // namespace csLie
