#include <csCore/testUtils/AssertCatch2.hpp>

#include <csLie/Types.hpp>
#include <csLie/SO2_impl.hpp>
#include <csLie/Constants.hpp>

#include <csLie/tests/TestUtils.hpp>

namespace csLie {
  

  TEMPLATE_TEST_CASE("SO2 test - base", "[SO2]", float32_t, float64_t) {
    using S = TestType;
    S const PRECISION = Constants<S>::dummy_precision;

    // create a identity rotation matrix
    SO2<S> const rI{ SO2<S>::Identity()};

    // create a zero algebra
    SO2Alg<S> const d{ SO2Alg<S>::Zero() };

    // call with Matrix<1,1>
    // exp the zero algebra, it has to be the identity matrix
    SO2<S> const exp0{ so2Exp(d) };
    REQUIRE((exp0 - rI).isZero());

    // TODO(sceriani) impl properly dynamic checks
    // so2Exp(Eigen::Vector2<S>()); // does not compile
    // so2Exp(Eigen::VectorX<S>()); // compile but crash

    // call with float
    // exp the 90deg (M_PI / 2), expect [0 -1; 1 0]
    SO2<S> const exp90{ so2Exp<S>(Constants<S>::pi / static_cast<S>(2.0)) };
    SO2<S> exp90_gt{};
    exp90_gt << 
      0, -1,
      1, 0;
    REQUIRE((exp90 - exp90_gt).isZero());


    Eigen::VectorX<S> var(5);
    var.setRandom();
    for (size_t i = 0; i < var.size(); i++) {
      auto retF = so2Exp(var(i));
      //auto retBlock = so2Exp(var.template segment<2>(i)); // does not compile
      auto retBlock = so2Exp(var.template segment<1>(i));
    }

    Eigen::MatrixX<S> m24(2, 4);
    so2Exp(static_cast<S>(0.), m24.template block<2, 2>(0, 0));
    so2Exp(0., m24.template block<2, 2>(0, 2));
    //so2Exp(0., m24.block<2, 1>(0, 2)); // does not compile

    so2Exp(var.template segment<1>(0), m24.template block<2, 2>(0, 0));
    // so2Exp(var.template segment<1>(1), m24.template block<2, 2>(0, 2));

    so2Exp(var.segment(0,1), m24.block(0, 0,2,2));

  }
}