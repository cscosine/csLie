project(csLieTest  LANGUAGES CXX)

cs_add_executable(csLieTest)

include(CTest)
cs_find_package(csLieTest PUBLIC Catch2 REQUIRED CONFIG)
include(Catch)

cs_target_sources(csLieTest 
	PRIVATE
	SE2.cpp
	S2.cpp
)

cs_target_link_libraries(csLieTest PRIVATE csLie::csLie Catch2::Catch2WithMain)
  
# add to tests
catch_discover_tests(csLieTest)

cs_install_target(csLieTest csLieTestExport DEFAULT)

cs_finalize_target(csLieTest)