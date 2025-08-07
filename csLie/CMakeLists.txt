project(csLie LANGUAGES CXX)


# prepare custom project asserts, note: use standard find_package
find_package(csCore REQUIRED COMPONENTS assert CONFIG)

list(APPEND CMAKE_MODULE_PATH "${csCore_BASE_DIR}/cmake/csCore/assert/")
include(csAssert)

# set defaults for Paranoid additional configuration
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT RELEASE Paranoid ON)
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT RELWITHDEBINFO Paranoid ON)
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT DEBUG Paranoid ON)
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT PARANOID Paranoid ON)
csSetAssertDefaultOnForConfig(${PROJECT_NAME} ASSERT SAFE_MATH Paranoid ON)

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
    	SE2.hpp
		SE2_impl.hpp
    	SE3.hpp
		SE3_impl.hpp
    	SO2.hpp
		SO2_impl.hpp
    	SO3.hpp
		SO3_impl.hpp
    	Utils.hpp
		Utils_impl.hpp
	PRIVATE
)

cs_target_link_libraries(csLie
	INTERFACE
		Eigen3::Eigen
    	csCore::core
		csCore::eigenUtils
		csLie::csLieAssert
		csLieSol::csLieSol_csCMakeBaseInclude
)

add_subdirectory(tests)

cs_install_target(csLie csLie DEFAULT)

cs_finalize_target(csLie)
