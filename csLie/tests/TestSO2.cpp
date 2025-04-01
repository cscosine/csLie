#include <csCore/testUtils/AssertCatch2.hpp>

#include <csLie/Types.hpp>
#include <csLie/SO2_impl.hpp>
#include <csLie/Constants.hpp>

#include <csLie/tests/TestUtils.hpp>

namespace csLie {
  

  TEMPLATE_TEST_CASE("SO2 test - base", "[SO2]", float32_t, float64_t) {
    using S = TestType;

    // create a identity rotation matrix
    SO2<S> const rI{ SO2<S>::Identity()};

    // create a zero algebra
    SO2Alg<S> const d0{ SO2Alg<S>::Zero() };

    // exp the zero algebra, it has to be the identity matrix
    SO2<S> const exp0{ so2Exp(d0) };
    REQUIRE((exp0 - rI).isZero());

    // exp the 90deg (M_PI / 2), expect [0 -1; 1 0]
    SO2<S> const exp90{ so2Exp<S>(Constants<S>::pi / static_cast<S>(2.0)) };
    SO2<S> exp90_gt{};
    exp90_gt << 
      0, -1,
      1, 0;
    REQUIRE((exp90 - exp90_gt).isZero());

  }

  TEMPLATE_TEST_CASE("SO2 test - variants", "[SO2]", float32_t, float64_t) {
    using S = TestType;

    /*
    so2Exp variants
    +------------+----------------------+
    | delta Type |  Return/Output Type  |
    |------------+----------------------|
    |   Matrix   |  return SO2          |
    |   Matrix   |  output matrix       |
    |   Matrix   |  output block        |
    |------------+----------------------|
    |   Scalar   |  return SO2          |
    |   Scalar   |  output matrix       |
    |   Scalar   |  output block        |
    |------------+----------------------|
    |   Block    |  return SO2          |
    |   Block    |  output matrix       |
    |   Block    |  output block        |
    +------------+----------------------+
    */

    SO2Alg<S> const deltaM{SO2Alg<S>::Random()};
    S const deltaS{deltaM(0)};
    SO2<S> retM = so2Exp(deltaM);
    SO2<S> retS = so2Exp(deltaS);

    SO2<S> outM;
    so2Exp(deltaM, outM);
    SO2<S> outS; 
    so2Exp(deltaS, outS);
    
    REQUIRE((outM - outS).isZero());

    // static blocks
    Eigen::Matrix<S,2,4> retB{};
    so2Exp(deltaM, retB.template block<2,2>(0,0));
    REQUIRE((retB.template block<2,2>(0,0) - outS).isZero());
    so2Exp(deltaS, retB.template block<2,2>(0,2));
    REQUIRE((retB.template block<2,2>(0,2) - outS).isZero());

    // dynamic blocks
    Eigen::Matrix<S,2,4> retBd{};
    so2Exp(deltaM, retBd.block(0,0,2,2));
    REQUIRE((retBd.template block<2,2>(0,0) - outS).isZero());
    so2Exp(deltaS, retBd.block(0,2,2,2));
    REQUIRE((retBd.template block<2,2>(0,2) - outS).isZero());

    // block in input (static and dynamic)
    static constexpr Eigen::Index N{5};
    Eigen::Vector<S,N> deltaV{Eigen::Vector<S,N>::Random()};
    std::array<SO2<S>, N> retV{};
    std::array<SO2<S>, N> retVd{};
    std::array<SO2<S>, N> outV{};
    std::array<SO2<S>, N> outVd{};

    Eigen::Matrix<S,2,N*2> outB{};
    Eigen::Matrix<S,2,N*2> outBd{};
    Eigen::Matrix<S,2,N*2> outBdd{};
    for(size_t i = 0; i < N; i++) {
      SO2<S> retI = so2Exp(deltaV(i)); 
      // return
      retV[i] = so2Exp(deltaV.template segment<1>(i));
      REQUIRE((retV[i] - retI).isZero());
      retVd[i] = so2Exp(deltaV.template segment(i,1));
      REQUIRE((retVd[i] - retI).isZero());

      // output
      so2Exp(deltaV.template segment<1>(i), outV[i]);
      REQUIRE((outV[i] - retI).isZero());
      so2Exp(deltaV.segment(i,1), outVd[i]);
      REQUIRE((outVd[i] - retI).isZero());

      // output block
      so2Exp(deltaV.template segment<1>(i), outB.template block<2,2>(0,2*i));
      REQUIRE((outB.template block<2,2>(0,2*i) - retI).isZero());
      so2Exp(deltaV.segment(i,1), outBd.template block<2,2>(0,2*i));
      REQUIRE((outBd.template block<2,2>(0,2*i) - retI).isZero());
      so2Exp(deltaV.segment(i,1), outBdd.block(0,2*i,2,2));
      REQUIRE((outBdd.block(0,2*i,2,2) - retI).isZero());
    }
  }

  TEMPLATE_TEST_CASE("SO2 test - delta type matrix, wrong dynamic sizes", "[SO2]", float32_t, float64_t) 
  {
    /*
    +------------+----------------------+
    | delta Type |  Return/Output Type  |
    |------------+----------------------|
    |   Matrix   |  return SO2          | test wrong delta matrix size
    |   Matrix   |  output matrix       | test wrong delta matrix size and wrong output matrix size
    |   Matrix   |  output block        | test wrong delta matrix size and wrong output block size
    +------------+----------------------+
    */
    using S = TestType;

    // so2Exp(Eigen::Matrix<S,3,1>()); // does not compile
    if(csLieAsserts::isAssertDebugEnabled()) {
      REQUIRE_ASSERT(so2Exp(Eigen::MatrixX<S>(1,2))); // throws, wrong delta
    }
    
    // Eigen::Matrix<S,3,3> ret33;
    // so2Exp(Eigen::Matrix<S,1,1>(), ret33); // does not compile
    if(csLieAsserts::isAssertDebugEnabled()) {
      Eigen::Matrix<S,2,2> ret22;
      REQUIRE_ASSERT(so2Exp(Eigen::MatrixX<S>(1,2), ret22)); // throws, wrong delta

      Eigen::MatrixX<S> ret33(3,3);
      REQUIRE_ASSERT(so2Exp(Eigen::MatrixX<S>(1,1), ret33)); // throws, wrong output
    }

    // so2Exp(Eigen::Matrix<S,1,1>(), ret24.template block<2,3>(0,0)); // does not compile
    if(csLieAsserts::isAssertDebugEnabled()) {
      Eigen::Matrix<S,2,4> ret24;
      REQUIRE_ASSERT(so2Exp(Eigen::MatrixX<S>(1,2), ret24.block(0,0, 2, 2))); // throws, wrong delta
      REQUIRE_ASSERT(so2Exp(Eigen::Matrix<S,1,1>(), ret24.block(0,0, 2, 3))); // throws, wrong output

    }
  }

  TEMPLATE_TEST_CASE("SO2 test - delta type scalar, wrong dynamic sizes", "[SO2]", float32_t, float64_t) 
  {
    /*
    +------------+----------------------+
    | delta Type |  Return/Output Type  |
    |------------+----------------------|  
    |   Scalar   |  return SO2          | no possible throw paths
    |   Scalar   |  output matrix       | test wrong output matrix size
    |   Scalar   |  output block        | test wrong output block size
    |------------+----------------------|
    */
    
    using S = TestType;

    if(csLieAsserts::isAssertDebugEnabled()) {
      Eigen::MatrixX<S> ret33(3,3);
      REQUIRE_ASSERT(so2Exp(S{0.0}, ret33)); // throws, wrong output
      REQUIRE_ASSERT(so2Exp(S{0.0}, ret33.block(0,0,2,3))); // throws, wrong output
    }
  }

  TEMPLATE_TEST_CASE("SO2 test - delta type block, wrong dynamic sizes", "[SO2]", float32_t, float64_t) 
  {
    /*
    +------------+----------------------+
    | delta Type |  Return/Output Type  |
    |------------+----------------------|  
    |   Block    |  return SO2          | wrong delta block size
    |   Block    |  output matrix       | wrong delta block size and output matrix size
    |   Block    |  output block        | wrong delta block size and output block size
    +------------+----------------------+
    */
    
    using S = TestType;

    if(csLieAsserts::isAssertDebugEnabled()) {
      Eigen::Matrix<S,2,5> input;
      REQUIRE_ASSERT(so2Exp(input.block(0,0,2,2))); // throws, wrong input block
      
      Eigen::MatrixX<S> output22(2,2);
      REQUIRE_ASSERT(so2Exp(input.block(0,0,2,2), output22)); // throws, wrong input block
      
      Eigen::MatrixX<S> output33(3,3);
      REQUIRE_ASSERT(so2Exp(input.block(0,0,1,1), output33)); // throws, wrong output

      REQUIRE_ASSERT(so2Exp(input.block(0,0,2,2), output33.block(0,0,2,2))); // throws, wrong input block
      REQUIRE_ASSERT(so2Exp(input.block(0,0,1,1), output33.block(0,0,2,3))); // throws, wrong output block
    }
  }

}