// Accuracy of the SO3 functions and of their coefficients against exact results (reference/data/SO3.hpp, computed
// in quad precision), in units of the epsilon of the tested precision. The inputs cover 0, tiny angles, 1 degree, the
// switch of the formulas, pi - 1e-k, exactly pi and beyond (to 2 pi - 1e-3 and 10), on axis-aligned and generic axes.
#include "reference/data/SO3.hpp"

#include <cmath>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SO3_impl.hpp>
#include <csLie/Types.hpp>

namespace csLie {

namespace {

/// the maximum errors accepted, in epsilon of the tested precision
constexpr double so3LogTolerance = 3.0;
// the exp: R = I + A W + B W^2 sums terms up to 2 in magnitude (B theta^2 <= 2), each a few roundings
constexpr double so3ExpTolerance = 4.0;
constexpr double coeffTolerance = 3.0;
// D = (B(x) - C(x)) / (4 A(x)) at x = theta / 2: beyond theta = 2 pi B(x) - C(x) cancels, by 2.5 at theta = 10
constexpr double coeffDTolerance = 4.0;
// the jacobians, J = A I + B W + C omega omega^T and J^-1 = K I - W / 2 + D omega omega^T: relative to the largest
// entry (J^-1 grows as 1 / (2 pi - theta) near 2 pi, where it is singular)
constexpr double so3JTolerance = 4.0;

template <class S> auto const &so3LogCases() {
  if constexpr (std::is_same_v<S, float32_t>) {
    return testReference::so3LogCasesF32;
  } else {
    return testReference::so3LogCasesF64;
  }
}

template <class S> auto const &so3ExpCases() {
  if constexpr (std::is_same_v<S, float32_t>) {
    return testReference::so3ExpCasesF32;
  } else {
    return testReference::so3ExpCasesF64;
  }
}

template <class S> auto const &so3JlCases() {
  if constexpr (std::is_same_v<S, float32_t>) {
    return testReference::so3JlCasesF32;
  } else {
    return testReference::so3JlCasesF64;
  }
}

template <class S> auto const &so3JlInvCases() {
  if constexpr (std::is_same_v<S, float32_t>) {
    return testReference::so3JlInvCasesF32;
  } else {
    return testReference::so3JlInvCasesF64;
  }
}

template <class S> auto const &coeffCases() {
  if constexpr (std::is_same_v<S, float32_t>) {
    return testReference::coeffCasesF32;
  } else {
    return testReference::coeffCasesF64;
  }
}

/// relative error of a coefficient, in epsilon of S
template <class S> double relError(S const got, double const exact) {
  return std::abs(static_cast<double>(got) - exact) / std::abs(exact) / static_cast<double>(Constants<S>::epsilon);
}

double norm3(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z); }

/// relative error of a log (absolute when the exact log is zero), in epsilon of S; near pi the log is defined up to
/// w -> w (1 - 2 pi / |w|), the closest of the two is taken
template <class S> double so3LogError(SO3Alg<S> const &got, double const (&exact)[3]) {
  double const g[3] = {static_cast<double>(got(0)), static_cast<double>(got(1)), static_cast<double>(got(2))};
  double const n = norm3(exact[0], exact[1], exact[2]);
  double err = norm3(g[0] - exact[0], g[1] - exact[1], g[2] - exact[2]) / (n > 0 ? n : 1.0);
  constexpr double pi = Constants<double>::pi;
  if (n > pi - 0.1) {
    double const f = 1.0 - 2.0 * pi / n;
    err = std::min(err, norm3(g[0] - f * exact[0], g[1] - f * exact[1], g[2] - f * exact[2]) / n);
  }
  return err / static_cast<double>(Constants<S>::epsilon);
}

/// the largest error of the entries of a 3x3 matrix, relative to its largest exact entry (to 1 when smaller), in
/// epsilon of S
template <class S> double matrixError(Eigen::Matrix<S, 3, 3> const &got, double const (&exact)[9]) {
  double err = 0;
  double scale = 1;
  for (int i = 0; i < 9; ++i) {
    err = std::max(err, std::abs(static_cast<double>(got(i / 3, i % 3)) - exact[i]));
    scale = std::max(scale, std::abs(exact[i]));
  }
  return err / scale / static_cast<double>(Constants<S>::epsilon);
}

/// the condition number of J^-1 in theta: near 2 pi an error of epsilon in theta (computed as |omega|, rounded)
/// moves J^-1 by theta / |2 pi - theta| epsilon, relative
double so3JinvConditioning(double const (&omega)[3]) {
  double const theta = norm3(omega[0], omega[1], omega[2]);
  return std::max(1.0, theta / std::abs(2.0 * Constants<double>::pi - theta));
}

/// checks fn (so3Exp, so3J, so3Jinv) on every case: the error, divided by the condition number of the case
/// (conditioning(omega)), within tolerance
template <class S, class Cases, class Fn, class Conditioning>
void checkMatrixCases(Cases const &cases, Fn fn, Conditioning conditioning, double const tolerance) {
  double maxErr = 0;
  for (auto const &c : cases) {
    SO3Alg<S> const w(static_cast<S>(c.omega[0]), static_cast<S>(c.omega[1]), static_cast<S>(c.omega[2]));
    Eigen::Matrix<S, 3, 3> const M = fn(w);
    double const err = matrixError(M, c.M) / conditioning(c.omega);
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": error " << err << " eps");
    CHECK(M.allFinite());
    CHECK(err <= tolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

} // namespace

TEMPLATE_TEST_CASE("SO3 log accuracy", "[SO3][accuracy]", float32_t, float64_t) {
  using S = TestType;
  double maxErr = 0;
  for (auto const &c : so3LogCases<S>()) {
    SO3<S> R;
    for (int i = 0; i < 9; ++i)
      R(i / 3, i % 3) = static_cast<S>(c.R[i]);
    SO3Alg<S> const w = so3Log(R);
    double const err = so3LogError(w, c.omega);
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": log = " << w.transpose() << ", error " << err << " eps");
    CHECK(w.allFinite());
    CHECK(err <= so3LogTolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

TEMPLATE_TEST_CASE("SO3 exp accuracy", "[SO3][accuracy]", float32_t, float64_t) {
  using S = TestType;
  // the entries are at most 1: the absolute error
  checkMatrixCases<S>(
      so3ExpCases<S>(), [](SO3Alg<S> const &w) { return so3Exp(w); }, [](double const(&)[3]) { return 1.0; },
      so3ExpTolerance);
}

TEMPLATE_TEST_CASE("SO3 left jacobian accuracy", "[SO3][accuracy]", float32_t, float64_t) {
  using S = TestType;
  checkMatrixCases<S>(
      so3JlCases<S>(), [](SO3Alg<S> const &w) { return so3J(w); }, [](double const(&)[3]) { return 1.0; },
      so3JTolerance);
}

TEMPLATE_TEST_CASE("SO3 inverse left jacobian accuracy", "[SO3][accuracy]", float32_t, float64_t) {
  using S = TestType;
  checkMatrixCases<S>(
      so3JlInvCases<S>(), [](SO3Alg<S> const &w) { return so3Jinv(w); }, so3JinvConditioning, so3JTolerance);
}

TEMPLATE_TEST_CASE("SO3 exp coefficients accuracy", "[SO3][accuracy]", float32_t, float64_t) {
  using S = TestType;
  for (auto const &c : coeffCases<S>()) {
    _private::ExpCoeffs<S> const k = _private::expCoeffs(static_cast<S>(c.theta));
    double const errA = relError(k.A, c.A);
    double const errB = relError(k.B, c.B);
    INFO(c.label << ": A error " << errA << " eps, B error " << errB << " eps");
    CHECK(errA <= coeffTolerance);
    CHECK(errB <= coeffTolerance);
  }
}

TEMPLATE_TEST_CASE("SO3 jacobian coefficients accuracy", "[SO3][accuracy]", float32_t, float64_t) {
  using S = TestType;
  for (auto const &c : coeffCases<S>()) {
    _private::JCoeffs<S> const k = _private::jCoeffs(static_cast<S>(c.theta));
    double const errA = relError(k.A, c.A);
    double const errB = relError(k.B, c.B);
    double const errC = relError(k.C, c.C);
    double const errD = relError(k.D, c.D);
    double const errK = relError(k.K, c.K);
    INFO(c.label << ": A " << errA << " eps, B " << errB << " eps, C " << errC << " eps, D " << errD << " eps, K "
                 << errK << " eps");
    CHECK(errA <= coeffTolerance);
    CHECK(errB <= coeffTolerance);
    CHECK(errC <= coeffTolerance);
    CHECK(errD <= coeffDTolerance);
    CHECK(errK <= coeffTolerance);
  }
}

} // namespace csLie
