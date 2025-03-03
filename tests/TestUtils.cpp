#include <tests/TestUtils.hpp>

#include <Eigen/Dense>
#include <csCore/BasicTypes.hpp>

#include <catch2/catch_all.hpp>

namespace csLie {

  template <class S> void checkValidRotationMatrix(Eigen::Matrix3<S> const& rot) {
    S const PRECISION = Eigen::NumTraits<S>::dummy_precision();
    // det = 1
    REQUIRE_THAT(rot.determinant(), Catch::Matchers::WithinAbs(1.0, PRECISION));
    // cols are orthogonal
    REQUIRE_THAT(rot.col(0).dot(rot.col(1)), Catch::Matchers::WithinAbs(0.0, PRECISION));
    REQUIRE_THAT(rot.col(0).dot(rot.col(2)), Catch::Matchers::WithinAbs(0.0, PRECISION));
    REQUIRE_THAT(rot.col(1).dot(rot.col(2)), Catch::Matchers::WithinAbs(0.0, PRECISION));
    // rows are orthogonal
    REQUIRE_THAT(rot.row(0).dot(rot.row(1)), Catch::Matchers::WithinAbs(0.0, PRECISION));
    REQUIRE_THAT(rot.row(0).dot(rot.row(2)), Catch::Matchers::WithinAbs(0.0, PRECISION));
    REQUIRE_THAT(rot.row(1).dot(rot.row(2)), Catch::Matchers::WithinAbs(0.0, PRECISION));
    // cols are unit vectors
    REQUIRE_THAT(rot.col(0).norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));
    REQUIRE_THAT(rot.col(1).norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));
    REQUIRE_THAT(rot.col(2).norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));
    // rows are unit vectors
    REQUIRE_THAT(rot.row(0).norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));
    REQUIRE_THAT(rot.row(1).norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));
    REQUIRE_THAT(rot.row(2).norm(), Catch::Matchers::WithinAbs(1.0, PRECISION));
    // x cross y = z
    REQUIRE((rot.col(0).cross(rot.col(1)) - rot.col(2)).isZero());
    // y cross z = x
    REQUIRE((rot.col(1).cross(rot.col(2)) - rot.col(0)).isZero());
    // z cross x = y
    REQUIRE((rot.col(2).cross(rot.col(0)) - rot.col(1)).isZero());
    // R^T * R = I
    REQUIRE((rot * rot.transpose() - Eigen::Matrix3<S>::Identity()).isZero());
    REQUIRE((rot.transpose() * rot - Eigen::Matrix3<S>::Identity()).isZero());
  }

  template void checkValidRotationMatrix<float32_t>(Eigen::Matrix3<float32_t> const& rot);
  template void checkValidRotationMatrix<float64_t>(Eigen::Matrix3<float64_t> const& rot);

}