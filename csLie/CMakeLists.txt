project(csLie LANGUAGES CXX)


# prepare custom project asserts, note: use standard find_package
find_package(csCore REQUIRED COMPONENTS assert CONFIG)

list(APPEND CMAKE_MODULE_PATH "${csCore_BASE_DIR}/cmake/csCore/assert/")
include(csAssert)
csAssertGenerateProject(csLie csLie) #generate target csLie::csLieAssert
list(POP_BACK CMAKE_MODULE_PATH)

# prepare standard library
cs_add_library(csLie INTERFACE)
cs_add_library(csLie::csLie ALIAS csLie)

cs_find_package(csLie PUBLIC Eigen3 REQUIRED CONFIG)
cs_find_package(csLie PUBLIC csCore COMPONENTS core eigenUtils REQUIRED CONFIG)

cs_target_sources(csLie 
	PRIVATE_INSTALL 
		Types.hpp
		Constants.hpp
    S2.hpp
		S2_impl.hpp
    SO2.hpp
		SO2_impl.hpp
	PRIVATE
)

cs_target_link_libraries(csLie
	INTERFACE
		Eigen3::Eigen
    	csCore::core
		csCore::eigenUtils
		csLie::csLieAssert
		csLieLib::csLieLib_csCMakeBaseInclude
)

add_subdirectory(tests)

cs_install_target(csLie csLieExport DEFAULT)

cs_finalize_target(csLie)