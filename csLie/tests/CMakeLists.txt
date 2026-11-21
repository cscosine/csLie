project(csLieTest LANGUAGES CXX)

cs_add_executable(csLieTest)

cs_find_package(csLieTest PRIVATE csCore COMPONENTS testUtils math REQUIRED CONFIG)

include(CTest)
cs_find_package(csLieTest PRIVATE Catch2 3 REQUIRED CONFIG)
include(Catch)

cs_target_sources(
  csLieTest
  PRIVATE TestSE2.cpp TestS2.cpp TestSO2.cpp TestSO3.cpp TestSE3.cpp TestUtils.cpp
)

cs_target_link_libraries(csLieTest PRIVATE csLie::csLie csCore::math csCore::testUtils Catch2::Catch2WithMain)

# add to tests
catch_discover_tests(csLieTest DISCOVERY_MODE "PRE_TEST")
