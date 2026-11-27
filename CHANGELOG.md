# Changelog

The changes of each csLie release that matter to the code using it. csLie is released with csfoundation.

## v0.2.0 (unreleased)

Accuracy: every exp, log and jacobian is accurate to a few epsilon at every angle, in `float32_t` and `float64_t`,
tested against results computed in quad precision (`csLie/tests/reference/`). Before, several lost most of their
digits at small angles (e.g. `so3J` and `se3Exp` up to 2e7 epsilon in double at 1e-8 rad, `se2Exp` alike) or near pi
(`so3Jinv`), and `so3Log` gave NaN at pi.

- **Changed, breaking: names.** What to change in the code using csLie:

  | Before | Now |
  |---|---|
  | `so3J`, `so3Jinv` | `so3LeftJacobian`, `so3LeftJacobianInverse` |
  | `so3ToRPY`, `rpyToSo3` | `so3ToRollPitchYaw`, `rollPitchYawToSo3` |
  | `skewSymm` | `hat` |
  | `xyYToSe2`, `se2ToxyY` | `xyYawToSe2`, `se2ToXyYaw` |
  | `s2oplus`, `s2_d_oplus_d_delta_0`, `s2toZAlignedSO3` | `s2Oplus`, `s2OplusJacobianAtZero`, `s2Frame` |
  | `se3JacobianLeftZeroPerturbationPointTransformation` | `se3PointJacobianLeft` |
  | `normalize3dVectorWithJacobians`, `UnitVector3dWithJacobians` (`vo`, `d_vo_d_v`) | `normalizeWithJacobian`, `NormalizedWithJacobian` (`unit`, `jacobian`) |
  | types `SO2Alg`, `SE2Alg`, `SO3Alg`, `SE3Alg`, `S2Alg` (and `d` / `f`) | `SO2Tangent`, `SE2Tangent`, `SO3Tangent`, `SE3Tangent`, `S2Tangent` |
  | type `SE3Adj` | `SE3Adjoint` |
  | types `SO3J`, `SO3Jinv`, `s2oplusJac` | removed: `Eigen::Matrix3<S>`, `Eigen::Matrix<S, 3, 2>` |

- **Changed, breaking: behaviour.**
  - `so3ToRollPitchYaw` at gimbal lock (pitch = +-pi / 2) puts the free angle in yaw (roll = 0) by default; before,
    in roll (yaw = 0): `so3ToRollPitchYaw<GimbalPolicy::kYawZero>` gives that. Away from it, the same angles, now
    accurate: `rollPitchYawToSo3` of the result gives the rotation back to a few epsilon at every pitch.
  - `s2Frame` (and so `s2Oplus`, `s2OplusJacobianAtZero`) near +-x (x^2 > 1/2) uses a second chart, the pivot axis z:
    a different tangent basis there, well conditioned everywhere.
- New: `so3RightJacobian`, `so3RightJacobianInverse`, `so3OplusLeft`, `so3OplusRight`, `se3OplusLeft`,
  `se3OplusRight`, `se3PointJacobianRight`, `vee`; `se2Log`; `xyYawToSe2(xy, yawRad)` and `xyYawToSe2(x, y, yawRad)`;
  `so3ToRollPitchYawChecked` (empty at gimbal lock) and `GimbalPolicy`; `s2OplusInverse` (the inverse of `s2Oplus`
  at a point).
- The functions read an Eigen expression argument once (an expression such as `Random()` gives new values at each
  read: `s2Exp` gave NaN).
- The library builds with strict warnings (warnings as errors); every public function is compiled in `float32_t` and
  `float64_t` under them by the tests.
