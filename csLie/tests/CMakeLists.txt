project(csLieTest LANGUAGES CXX)

# the tests: every warning (cscmake_warnings), not as errors; the library code is checked under strict warnings by
# csLieInstantiationCheck below
set(CSCMAKE_DEFAULT_COMPILE_PROFILES cscmake_default cscmake_warnings)

# every public function of csLie in float32_t and float64_t, under strict warnings (compiled, not run)
cs_add_library(csLieInstantiationCheck OBJECT)
cs_target_sources(csLieInstantiationCheck PRIVATE InstantiationCheck.cpp)
cs_target_link_libraries(csLieInstantiationCheck PRIVATE csLie::csLie)
cs_target_use_compile_profiles(csLieInstantiationCheck cscmake_strict_warnings)

cs_add_executable(csLieTest)

cs_find_package(csLieTest PRIVATE csCore COMPONENTS testUtils math REQUIRED CONFIG)

include(CTest)
cs_find_package(csLieTest PRIVATE Catch2 3 REQUIRED CONFIG)
include(Catch)

cs_target_sources(
  csLieTest
  PRIVATE
    TestSE2.cpp
    TestSE2Accuracy.cpp
    TestS2.cpp
    TestS2Accuracy.cpp
    TestSO2.cpp
    TestSO3.cpp
    TestSO3Accuracy.cpp
    TestSE3.cpp
    TestSE3Accuracy.cpp
    TestPerturbations.cpp
    TestUtils.cpp
)

cs_target_link_libraries(csLieTest PRIVATE csLie::csLie csCore::math csCore::testUtils Catch2::Catch2WithMain)

# add to tests
catch_discover_tests(csLieTest DISCOVERY_MODE "PRE_TEST")

# the generator of the reference data (reference/generateReferenceData.cpp), where __float128 and libquadmath exist
# (GCC on x86-64 Linux): a test checks that the committed reference/data/ folder is its output, a target regenerates it
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES quadmath)
check_cxx_source_compiles(
  [=[
  #include <quadmath.h>
  int main() { __float128 x = 1; return atanq(x) > 0 ? 0 : 1; }
  ]=]
  CSLIE_HAS_QUADMATH
)
unset(CMAKE_REQUIRED_LIBRARIES)

if(CSLIE_HAS_QUADMATH)
  cs_add_executable(csLieGenerateReferenceData)
  cs_target_sources(csLieGenerateReferenceData PRIVATE reference/generateReferenceData.cpp)
  cs_target_link_libraries(csLieGenerateReferenceData PRIVATE quadmath)
  cs_add_executable(csLieCompareReferenceData)
  cs_target_sources(csLieCompareReferenceData PRIVATE reference/compareReferenceData.cpp)

  set(
    reference_args
    -DGENERATOR=$<TARGET_FILE:csLieGenerateReferenceData>
    -DCOMPARATOR=$<TARGET_FILE:csLieCompareReferenceData>
    -DREFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/reference/data
    -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/referenceData
  )
  add_test(
    NAME "csLie reference data up to date"
    COMMAND ${CMAKE_COMMAND} ${reference_args} -DMODE=CHECK -P ${CMAKE_CURRENT_SOURCE_DIR}/reference/ReferenceData.cmake
  )
  add_custom_target(
    csLieRegenerateReferenceData
    COMMAND ${CMAKE_COMMAND} ${reference_args} -DMODE=WRITE -P ${CMAKE_CURRENT_SOURCE_DIR}/reference/ReferenceData.cmake
    DEPENDS csLieGenerateReferenceData csLieCompareReferenceData
    COMMENT "Regenerating csLie/tests/reference/data/"
    VERBATIM
  )
else()
  message(STATUS "csLie: no __float128 / libquadmath, the reference data check is skipped")
endif()
