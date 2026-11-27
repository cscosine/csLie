// Accuracy of the SE3 exp and log against exact results (reference/data/SE3.hpp, computed in quad precision), in
// units of the epsilon of the tested precision. The inputs cover 0, tiny angles, 1 degree, pi - 1e-k, pi and beyond
// (the exp, to 10), axis-aligned and generic axes, unit to far translations.
#include "reference/data/SE3.hpp"

#include <cmath>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SE3_impl.hpp>
#include <csLie/Types.hpp>

namespace csLie {

namespace {

/// the maximum errors accepted, in epsilon of the tested precision
// the rotation entries; the translation relative to |rho| (|t| <= |rho|), to 1 when smaller
constexpr double se3ExpTolerance = 4.0;
// relative to |xi|, to 1 when smaller
constexpr double se3LogTolerance = 4.0;

template <class S> auto const &se3ExpCases() {
  if constexpr (std::is_same_v<S, float32_t>) {
    return testReference::se3ExpCasesF32;
  } else {
    return testReference::se3ExpCasesF64;
  }
}

template <class S> auto const &se3LogCases() {
  if constexpr (std::is_same_v<S, float32_t>) {
    return testReference::se3LogCasesF32;
  } else {
    return testReference::se3LogCasesF64;
  }
}

double norm(double const *v, int n) {
  double s = 0;
  for (int i = 0; i < n; ++i)
    s += v[i] * v[i];
  return std::sqrt(s);
}

} // namespace

TEMPLATE_TEST_CASE("SE3 exp accuracy", "[SE3][accuracy]", float32_t, float64_t) {
  using S = TestType;
  double maxErr = 0;
  for (auto const &c : se3ExpCases<S>()) {
    SE3Tangent<S> xi;
    for (int i = 0; i < 6; ++i)
      xi(i) = static_cast<S>(c.xi[i]);
    SE3<S> const T = se3Exp(xi);
    double const scale = std::max(1.0, norm(c.xi, 3));
    double err = 0;
    for (int i = 0; i < 16; ++i) {
      double const e = std::abs(static_cast<double>(T(i / 4, i % 4)) - c.T[i]);
      err = std::max(err, i % 4 == 3 ? e / scale : e);
    }
    err /= static_cast<double>(Constants<S>::epsilon);
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": error " << err << " eps");
    CHECK(T.matrix().allFinite());
    CHECK(err <= se3ExpTolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

TEMPLATE_TEST_CASE("SE3 log accuracy", "[SE3][accuracy]", float32_t, float64_t) {
  using S = TestType;
  double maxErr = 0;
  for (auto const &c : se3LogCases<S>()) {
    SE3Matrix<S> T;
    for (int i = 0; i < 16; ++i)
      T(i / 4, i % 4) = static_cast<S>(c.T[i]);
    SE3Tangent<S> const xi = se3Log(T);
    double d[6];
    for (int i = 0; i < 6; ++i)
      d[i] = static_cast<double>(xi(i)) - c.xi[i];
    double const err = norm(d, 6) / std::max(1.0, norm(c.xi, 6)) / static_cast<double>(Constants<S>::epsilon);
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": log = " << xi.transpose() << ", error " << err << " eps");
    CHECK(xi.allFinite());
    CHECK(err <= se3LogTolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

} // namespace csLie
