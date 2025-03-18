project(csLieTest  LANGUAGES CXX)

cs_add_executable(csLieTest)

include(CTest)
cs_find_package(csLieTest PUBLIC Catch2 3 REQUIRED CONFIG)
include(Catch)

cs_target_sources(csLieTest 
	PRIVATE
	# SE2.cpp
  TestUtils.hpp
  TestUtils.cpp
	TestS2.cpp
  TestSO2.cpp
)

cs_target_link_libraries(csLieTest 
  PRIVATE 
  csLie::csLie 
  csLieLib::csLieLib_csCMakeBaseInclude
  csLieLib::csLieLib_csCMakeCompileOptionsDefault
  Catch2::Catch2WithMain)
  
# add to tests
catch_discover_tests(csLieTest)

cs_install_target(csLieTest csLieTestExport)

cs_finalize_target(csLieTest)