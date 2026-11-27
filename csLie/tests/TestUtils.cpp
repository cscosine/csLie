#include <cmath>
#include <csCore/testUtils/AssertCatch2.hpp>
#include <csLie/Constants.hpp>
#include <csLie/Types.hpp>
#include <csLie/Utils_impl.hpp>

namespace csLie {

TEMPLATE_TEST_CASE("Utils test - normalize", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  Eigen::Vector3<S> const u(std::sqrt(S(3)) / S(3), std::sqrt(S(3)) / S(3), std::sqrt(S(3)) / S(3));
  UnitVector3dWithJacobians<S> ret{normalize3dVectorWithJacobians(u)};
  REQUIRE((ret.vo - u).isZero());
  Eigen::Matrix3<S> d_vo_d_v_gt;
  d_vo_d_v_gt << S(2) / S(3), S(-1) / S(3), S(-1) / S(3), S(-1) / S(3), S(2) / S(3), S(-1) / S(3), S(-1) / S(3),
      S(-1) / S(3), S(2) / S(3);
  REQUIRE((ret.d_vo_d_v - d_vo_d_v_gt).isZero());

  ret = normalize3dVectorWithJacobians(Eigen::Vector3<S>(2, 0, 0));
  REQUIRE((ret.vo - Eigen::Vector3<S>::UnitX()).isZero());
  d_vo_d_v_gt << 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5;
  REQUIRE((ret.d_vo_d_v - d_vo_d_v_gt).isZero());
}

TEMPLATE_TEST_CASE("Utils test - normalize, vector wrong size", "[SO3]", float32_t, float64_t) {
  using S = TestType;

  if (csLie::Asserts::isAssertDebugEnabled()) {
    REQUIRE_NO_ASSERT(normalize3dVectorWithJacobians(Eigen::Matrix<S, 3, 1>::UnitX()));
    REQUIRE_ASSERT(normalize3dVectorWithJacobians(Eigen::MatrixX<S>(1, 2))); // throws, wrong size
    REQUIRE_ASSERT(normalize3dVectorWithJacobians(Eigen::MatrixX<S>(2, 1))); // throws, wrong size

    Eigen::VectorX<S> v10 = Eigen::VectorX<S>::Zero(10);
    v10(0) = 1;

    REQUIRE_NO_ASSERT(normalize3dVectorWithJacobians(Eigen::MatrixX<S>::Identity(5, 7).block(0, 0, 3, 1)));
    REQUIRE_NO_ASSERT(normalize3dVectorWithJacobians(v10.segment(0, 3)));
    REQUIRE_ASSERT(normalize3dVectorWithJacobians(Eigen::MatrixX<S>::Identity(2, 2).block(0, 0, 1, 2)));
    REQUIRE_ASSERT(normalize3dVectorWithJacobians(Eigen::MatrixX<S>::Identity(2, 2).block(0, 0, 2, 1)));
    REQUIRE_ASSERT(normalize3dVectorWithJacobians(v10.segment(0, 2)));
  }

  if (csLie::Asserts::isAssertParanoidEnabled()) {
    REQUIRE_ASSERT(normalize3dVectorWithJacobians(Eigen::Vector3<S>::Zero())); // throws, zero vector cannot be
                                                                               // normalized
  }
}

} // namespace csLie
