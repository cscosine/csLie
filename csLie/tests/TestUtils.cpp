#include <cmath>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/Types.hpp>
#include <csLie/Utils_impl.hpp>

namespace csLie {

TEMPLATE_TEST_CASE("Utils test - normalize", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  Eigen::Vector3<S> const u(std::sqrt(S(3)) / S(3), std::sqrt(S(3)) / S(3), std::sqrt(S(3)) / S(3));
  NormalizedWithJacobian<S> ret{normalizeWithJacobian(u)};
  REQUIRE((ret.unit - u).isZero());
  Eigen::Matrix3<S> jacobian_gt;
  jacobian_gt << S(2) / S(3), S(-1) / S(3), S(-1) / S(3), S(-1) / S(3), S(2) / S(3), S(-1) / S(3), S(-1) / S(3),
      S(-1) / S(3), S(2) / S(3);
  REQUIRE((ret.jacobian - jacobian_gt).isZero());

  ret = normalizeWithJacobian(Eigen::Vector3<S>(2, 0, 0));
  REQUIRE((ret.unit - Eigen::Vector3<S>::UnitX()).isZero());
  jacobian_gt << 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5;
  REQUIRE((ret.jacobian - jacobian_gt).isZero());
}

TEMPLATE_TEST_CASE("Utils test - normalize, vector wrong size", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(normalizeWithJacobian(Eigen::Matrix<S, 3, 1>::UnitX()));
    REQUIRE_ASSERT(normalizeWithJacobian(Eigen::MatrixX<S>(1, 2))); // throws, wrong size
    REQUIRE_ASSERT(normalizeWithJacobian(Eigen::MatrixX<S>(2, 1))); // throws, wrong size

    Eigen::VectorX<S> v10 = Eigen::VectorX<S>::Zero(10);
    v10(0) = 1;

    REQUIRE_NO_ASSERT(normalizeWithJacobian(Eigen::MatrixX<S>::Identity(5, 7).block(0, 0, 3, 1)));
    REQUIRE_NO_ASSERT(normalizeWithJacobian(v10.segment(0, 3)));
    REQUIRE_ASSERT(normalizeWithJacobian(Eigen::MatrixX<S>::Identity(2, 2).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(normalizeWithJacobian(Eigen::MatrixX<S>::Identity(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(normalizeWithJacobian(v10.segment(0, 2)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    REQUIRE_ASSERT(normalizeWithJacobian(Eigen::Vector3<S>::Zero())); // throws, zero vector cannot be
                                                                      // normalized
  }
}

} // namespace csLie
