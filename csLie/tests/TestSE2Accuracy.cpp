// Accuracy of the SE2 exp and log against exact results (reference/data/SE2.hpp, computed in quad precision), in
// units of the epsilon of the tested precision. The inputs cover 0, tiny angles, 1 degree, pi - 1e-k, pi and beyond
// (the exp, to 10), both signs, unit to far translations.
#include "reference/data/SE2.hpp"

#include <cmath>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/SE2_impl.hpp>
#include <csLie/Types.hpp>

namespace csLie {

namespace {

/// the maximum errors accepted, in epsilon of the tested precision
constexpr double se2ExpTolerance = 3.0;
constexpr double se2LogTolerance = 3.0;

template <class S> auto const &se2ExpCases() {
  if constexpr (std::is_same_v<S, float32_t>) {
    return testReference::se2ExpCasesF32;
  } else {
    return testReference::se2ExpCasesF64;
  }
}

template <class S> auto const &se2LogCases() {
  if constexpr (std::is_same_v<S, float32_t>) {
    return testReference::se2LogCasesF32;
  } else {
    return testReference::se2LogCasesF64;
  }
}

double norm3(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z); }

} // namespace

TEMPLATE_TEST_CASE("SE2 exp accuracy", "[SE2][accuracy]", float32_t, float64_t) {
  using S = TestType;
  double maxErr = 0;
  for (auto const &c : se2ExpCases<S>()) {
    SE2Tangent<S> const xi(static_cast<S>(c.xi[0]), static_cast<S>(c.xi[1]), static_cast<S>(c.xi[2]));
    SE2<S> const T = se2Exp(xi);
    // the rotation entries are at most 1; the translation is relative to |rho| (|t| <= |rho|), to 1 when smaller
    double const scale = std::max(1.0, std::hypot(c.xi[0], c.xi[1]));
    double err = 0;
    for (int i = 0; i < 9; ++i) {
      double const e = std::abs(static_cast<double>(T(i / 3, i % 3)) - c.T[i]);
      err = std::max(err, i % 3 == 2 ? e / scale : e);
    }
    err /= static_cast<double>(Constants<S>::epsilon);
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": error " << err << " eps");
    CHECK(T.matrix().allFinite());
    CHECK(err <= se2ExpTolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

TEMPLATE_TEST_CASE("SE2 log accuracy", "[SE2][accuracy]", float32_t, float64_t) {
  using S = TestType;
  double maxErr = 0;
  for (auto const &c : se2LogCases<S>()) {
    SE2Matrix<S> T;
    for (int i = 0; i < 9; ++i)
      T(i / 3, i % 3) = static_cast<S>(c.T[i]);
    SE2Tangent<S> const xi = se2Log(T);
    // relative to |xi|, to 1 when smaller
    double const scale = std::max(1.0, norm3(c.xi[0], c.xi[1], c.xi[2]));
    double const err = norm3(static_cast<double>(xi(0)) - c.xi[0], static_cast<double>(xi(1)) - c.xi[1],
                             static_cast<double>(xi(2)) - c.xi[2]) /
                       scale / static_cast<double>(Constants<S>::epsilon);
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": log = " << xi.transpose() << ", error " << err << " eps");
    CHECK(xi.allFinite());
    CHECK(err <= se2LogTolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

} // namespace csLie
