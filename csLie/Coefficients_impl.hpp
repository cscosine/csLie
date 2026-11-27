// The coefficients of the exp, log and jacobians of SO3, SE2 and SE3 as functions of the rotation angle, accurate to a
// few epsilon at every angle: shared by SO3_impl.hpp, SE2_impl.hpp and SE3_impl.hpp.
#pragma once

#include <cmath>
#include <csLie/Constants.hpp>

namespace csLie {

namespace _private {

/// @brief the coefficients of the exp at the angle theta: SO3, R = I + A W + B W^2 (W = hat(omega),
/// theta = |omega|); SE2, t = [A, -B theta; B theta, A] rho
template <class S> struct ExpCoeffs {
  S A; ///< sin(theta) / theta
  S B; ///< (1 - cos(theta)) / theta^2
};

/// @brief the coefficients of the exp at the angle theta >= 0, from the half angle: with
/// h = sin(theta / 2) / (theta / 2) and c = cos(theta / 2), A = h c and B = h^2 / 2. Accurate to a few epsilon at every
/// angle, with no series and no switch: sin(x) / x has no cancellation for any x > 0 (B as (1 - cos(theta)) / theta^2
/// loses all its digits below theta ~ sqrt(epsilon)); at theta = 0, h = 1.
template <class S> ExpCoeffs<S> expCoeffs(S const theta) {
  using std::cos;
  using std::sin;
  S const half = theta / S(2);
  S const h = half > S(0) ? sin(half) / half : S(1);
  return {h * cos(half), S(0.5) * h * h};
}

/// @brief the coefficients of the SO3 left jacobian and of its inverse, also of SE2 and SE3
/// (W = hat(omega), theta = |omega|):
/// J = I + B W + C W^2 = A I + B W + C omega omega^T, J^-1 = I - W / 2 + D W^2 = K I - W / 2 + D omega omega^T
template <class S> struct JCoeffs {
  S A; ///< sin(theta) / theta
  S B; ///< (1 - cos(theta)) / theta^2
  S C; ///< (theta - sin(theta)) / theta^3
  S D; ///< (1 - A / (2 B)) / theta^2
  S K; ///< A / (2 B) = (theta / 2) cot(theta / 2) = 1 - D theta^2
};

/// @brief the number of terms of the series of C(x) = (x - sin(x)) / x^3 = sum_k (-x^2)^k / (2k + 3)! that makes its
/// truncation error below epsilon / 8 of C at x <= 2 (C(2) > 1/8): the first omitted term 4^n / (2n + 3)! below
/// epsilon / 64. A bound, computed in S at compile time
template <class S> constexpr int cSeriesTerms() {
  S term = S(1) / S(6); // the first omitted term, 4^n / (2n + 3)!, at n = 0
  int n = 0;
  while (term > Constants<S>::epsilon / S(64)) {
    ++n;
    term *= S(4) / S((2 * n + 2) * (2 * n + 3));
  }
  return n;
}

/// @brief C(x) = (x - sin(x)) / x^3 at x >= 0, given sin(x): accurate to a few epsilon. The closed form loses digits
/// to the cancellation of x - sin(x) below x ~ 2: there the series, with as many terms as the precision of S needs
template <class S> S coeffC(S const x, S const sinX) {
  if (x > S(2)) {
    return (x - sinX) / (x * x * x);
  }
  constexpr int terms = cSeriesTerms<S>();
  // Horner: C = (1 - x^2 / (4 5) (1 - x^2 / (6 7) (1 - ...))) / 6
  S const x2 = x * x;
  S r = S(1);
  for (int k = terms - 1; k >= 1; --k) {
    r = S(1) - x2 * r / S((2 * k + 2) * (2 * k + 3));
  }
  return r / S(6);
}

/// @brief the coefficients of the left jacobian and of its inverse at the angle theta >= 0, accurate to a few
/// epsilon at every angle (J^-1 is singular at theta = 2 pi). From the half angle x = theta / 2, with no cancellation:
/// A and B as expCoeffs, C(theta) = (C(x) + A(x) B(x)) / 4, D(theta) = (B(x) - C(x)) / (4 A(x)), K = cos(x) / A(x)
/// (C(x) ~ 1/6 against A(x) B(x) ~ 1/2 and B(x) ~ 1/2; beyond theta = 2 pi, B(x) - C(x) cancels mildly, by 2.5 at
/// theta = 10). C and D in closed form, (theta - sin(theta)) / theta^3 and
/// (1 - A / (2 B)) / theta^2, lose all their digits at small angles.
template <class S> JCoeffs<S> jCoeffs(S const theta) {
  using std::cos;
  using std::sin;
  S const half = theta / S(2);
  S const quarter = theta / S(4);
  S const sinHalf = sin(half);
  S const cosHalf = cos(half);
  S const h = half > S(0) ? sinHalf / half : S(1);             // A(theta / 2)
  S const hq = quarter > S(0) ? sin(quarter) / quarter : S(1); // A(theta / 4)
  S const bHalf = S(0.5) * hq * hq;                            // B(theta / 2)
  S const cHalf = coeffC(half, sinHalf);                       // C(theta / 2)
  return {h * cosHalf, S(0.5) * h * h, (cHalf + h * bHalf) / S(4), (bHalf - cHalf) / (S(4) * h), cosHalf / h};
}

} // namespace _private

} // namespace csLie
