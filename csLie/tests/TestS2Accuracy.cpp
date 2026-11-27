// Accuracy of the S2 functions against exact results (reference/data/S2.hpp, computed in quad precision), in units of
// the epsilon of the tested precision. The inputs cover the axes, generic vectors, vectors near +-x (the second chart)
// and near the switch of the charts, angles from 0 to pi and beyond, and w = -v.
#include "reference/data/S2.hpp"

#include <cmath>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/S2_impl.hpp>
#include <csLie/Types.hpp>
#include <cstddef>

namespace csLie {

namespace {

/// the maximum errors accepted, in epsilon of the tested precision
constexpr double s2ExpTolerance = 3.0;
constexpr double s2FrameTolerance = 3.0;
constexpr double s2OplusTolerance = 4.0;
// relative to |delta|, times the condition number of the direction of delta, max(1, theta / sin(theta)) (the cut
// locus at pi)
constexpr double s2OplusInverseTolerance = 4.0;
// s2Oplus(v, s2OplusInverse(v, w)) against w
constexpr double s2RoundTripTolerance = 6.0;

#define CSLIE_CASES(fn, name)                                                                                          \
  template <class S> auto const &fn() {                                                                                \
    if constexpr (std::is_same_v<S, float32_t>) {                                                                      \
      return testReference::name##F32;                                                                                 \
    } else {                                                                                                           \
      return testReference::name##F64;                                                                                 \
    }                                                                                                                  \
  }
CSLIE_CASES(s2ExpCases, s2ExpCases)
CSLIE_CASES(s2FrameCases, s2FrameCases)
CSLIE_CASES(s2OplusCases, s2OplusCases)
CSLIE_CASES(s2OplusInverseCases, s2OplusInverseCases)
#undef CSLIE_CASES

template <class S, std::size_t N> Eigen::Matrix<S, static_cast<int>(N), 1> vec(double const (&a)[N]) {
  Eigen::Matrix<S, static_cast<int>(N), 1> o;
  for (std::size_t i = 0; i < N; ++i)
    o(static_cast<Eigen::Index>(i)) = static_cast<S>(a[i]);
  return o;
}

/// the largest absolute error of the entries, in epsilon of S
template <class S, int R, int C> double absError(Eigen::Matrix<S, R, C> const &got, double const *exact) {
  double err = 0;
  for (int i = 0; i < R * C; ++i)
    err = std::max(err, std::abs(static_cast<double>(got(i / C, i % C)) - exact[i]));
  return err / static_cast<double>(Constants<S>::epsilon);
}

} // namespace

TEMPLATE_TEST_CASE("S2 exp accuracy", "[S2][accuracy]", float32_t, float64_t) {
  using S = TestType;
  double maxErr = 0;
  for (auto const &c : s2ExpCases<S>()) {
    S2<S> const v = s2Exp(vec<S>(c.delta));
    double const err = absError(v, c.v);
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": error " << err << " eps");
    CHECK(err <= s2ExpTolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

TEMPLATE_TEST_CASE("S2 frame accuracy", "[S2][accuracy]", float32_t, float64_t) {
  using S = TestType;
  double maxErr = 0;
  for (auto const &c : s2FrameCases<S>()) {
    SO3<S> const R = s2Frame(vec<S>(c.v));
    double const err = absError(R, c.R); // row major in the table
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": error " << err << " eps");
    CHECK(R.allFinite());
    CHECK(err <= s2FrameTolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

TEMPLATE_TEST_CASE("S2 oplus accuracy", "[S2][accuracy]", float32_t, float64_t) {
  using S = TestType;
  double maxErr = 0;
  for (auto const &c : s2OplusCases<S>()) {
    S2<S> const w = s2Oplus(vec<S>(c.v), vec<S>(c.delta));
    double const err = absError(w, c.w);
    maxErr = std::max(maxErr, err);
    INFO(c.label << ": error " << err << " eps");
    CHECK(err <= s2OplusTolerance);
  }
  INFO("max error " << maxErr << " eps");
  SUCCEED();
}

TEMPLATE_TEST_CASE("S2 oplus inverse accuracy", "[S2][accuracy]", float32_t, float64_t) {
  using S = TestType;
  double maxErr = 0;
  double maxRoundTrip = 0;
  for (auto const &c : s2OplusInverseCases<S>()) {
    S2<S> const v = vec<S>(c.v);
    S2<S> const w = vec<S>(c.w);
    S2Tangent<S> const delta = s2OplusInverse(v, w);
    double const theta = std::hypot(c.delta[0], c.delta[1]);
    double const conditioning = theta > 0 ? std::max(1.0, theta / std::sin(theta)) : 1.0;
    double const err =
        std::hypot(static_cast<double>(delta(0)) - c.delta[0], static_cast<double>(delta(1)) - c.delta[1]) /
        std::max(theta, 1e-300) / conditioning / static_cast<double>(Constants<S>::epsilon);
    bool const antipodal = std::fpclassify(c.w[0] + c.v[0]) == FP_ZERO && std::fpclassify(c.w[1] + c.v[1]) == FP_ZERO &&
                           std::fpclassify(c.w[2] + c.v[2]) == FP_ZERO;
    double const roundTrip = absError(S2<S>(s2Oplus(v, delta)), c.w);
    maxErr = std::max(maxErr, antipodal ? 0.0 : err);
    maxRoundTrip = std::max(maxRoundTrip, roundTrip);
    INFO(c.label << ": delta " << delta.transpose() << ", error " << err << " eps, round trip " << roundTrip << " eps");
    CHECK(delta.allFinite());
    CHECK(roundTrip <= s2RoundTripTolerance);
    if (antipodal) {
      // the documented direction, b1
      CHECK(std::fpclassify(delta(0) - Constants<S>::pi) == FP_ZERO);
      CHECK(std::fpclassify(delta(1)) == FP_ZERO);
    } else if (theta > 0) {
      CHECK(err <= s2OplusInverseTolerance);
    } else {
      CHECK(delta.isZero(0));
    }
  }
  INFO("max error " << maxErr << " eps, max round trip " << maxRoundTrip << " eps");
  SUCCEED();
}

} // namespace csLie
