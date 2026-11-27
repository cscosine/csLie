project(csLie LANGUAGES CXX)

# prepare custom project asserts, note: use standard find_package
find_package(csCore REQUIRED COMPONENTS assert CONFIG)

list(APPEND CMAKE_MODULE_PATH "${csCore_BASE_DIR}/cmake/csCore/assert/")
include(csAssert)

# set defaults for Paranoid additional configuration
cs_assert_set_default_on_for_config(${PROJECT_NAME} ASSERT RELEASE Paranoid ON)
cs_assert_set_default_on_for_config(${PROJECT_NAME} ASSERT RELWITHDEBINFO Paranoid ON)
cs_assert_set_default_on_for_config(${PROJECT_NAME} ASSERT DEBUG Paranoid ON)
cs_assert_set_default_on_for_config(${PROJECT_NAME} ASSERT PARANOID Paranoid ON)
cs_assert_set_default_on_for_config(${PROJECT_NAME} ASSERT SAFE_MATH Paranoid ON)

cs_assert_generate_project(csLie csLie) #generate target csLie::csLieAssert
list(POP_BACK CMAKE_MODULE_PATH)

# prepare standard library
cs_add_library(csLie INTERFACE)
cs_add_library(csLie::csLie ALIAS csLie)

cs_find_package(csLie PUBLIC Eigen3 REQUIRED CONFIG)
cs_find_package(csLie PUBLIC csCore COMPONENTS core eigenUtils REQUIRED CONFIG)

cs_target_sources(
  csLie
  PRIVATE_INSTALL
    Types.hpp
    Coefficients_impl.hpp
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
)

cs_target_link_libraries(
  csLie
  INTERFACE Eigen3::Eigen csCore::core csCore::eigenUtils csLie::csLieAssert csLie::cscmake_base_include
)

add_subdirectory(tests)

cs_install_target(csLie COMPONENT csLie DEFAULT)

# compile each single header in standalone, ensures includes are self contained
cs_target_check_installed_headers(csLie COMPILE_PROFILES cscmake_strict_warnings csLie_paranoid)
