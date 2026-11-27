// Accuracy of the SO3 functions and of their coefficients against exact results (reference/data/SO3.hpp, computed
// in quad precision), in units of the epsilon of the tested precision. The inputs cover 0, tiny angles, 1 degree, the
// switch of the formulas, pi - 1e-k and exactly pi, on axis-aligned and generic axes.
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
  double maxErr = 0;
  for (auto const &c : so3ExpCases<S>()) {
    SO3Alg<S> const w(static_cast<S>(c.omega[0]), static_cast<S>(c.omega[1]), static_cast<S>(c.omega[2]));
    SO3<S> const R = so3Exp(w);
    // the entries are at most 1: the absolute error, in epsilon of S
    double err = 0;
    for (int i = 0; i < 9; ++i)
      err = std::max(err, std::abs(static_cast<double>(R(i / 3, i % 3)) - c.R[i]));
    err /= static_cast<double>(Constants<S>::epsilon);
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": error " << err << " eps");
    CHECK(R.allFinite());
    CHECK(err <= so3ExpTolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

TEMPLATE_TEST_CASE("SO3 exp coefficients accuracy", "[SO3][accuracy]", float32_t, float64_t) {
  using S = TestType;
  for (auto const &c : coeffCases<S>()) {
    _private::SO3ExpCoeffs<S> const k = _private::so3ExpCoeffs(static_cast<S>(c.theta));
    double const errA = relError(k.A, c.A);
    double const errB = relError(k.B, c.B);
    INFO(c.label << ": A error " << errA << " eps, B error " << errB << " eps");
    CHECK(errA <= coeffTolerance);
    CHECK(errB <= coeffTolerance);
  }
}

} // namespace csLie
