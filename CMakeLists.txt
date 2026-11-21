cmake_minimum_required(VERSION 3.31)

# before project(): the oldest macOS the binaries run on (csCore)
include(csCMakeSetup)
cs_macos_deployment_target(13.3)

project(csLieSol VERSION 0.1.0 LANGUAGES CXX)

# add the csCMake system
include(csCMake)

# Paranoid: undefined behavior sanitizer, used by the installed headers check
cs_add_compile_profile_rule(PROFILE csLie_paranoid CONFIGS Paranoid PROFILES cscmake_sanitize_undefined)

cs_initialize_project(
  NAMESPACE csLie
  CXX_STANDARD 17
  # the default configurations supported by this project
  DEFAULT_CONFIG Debug Release
  DEFAULT_CONFIG_MAPPED RelWithDebInfo
  MAPPED_TO "RelWithDebInfo,Release,"
  # a custom configuration, based on Debug (imported targets: Debug, else Release)
  CUSTOM_CONFIG Paranoid
  BASED_ON "Debug,Release,"
  POSTFIX "dp"
)

# add libs, executables, etc
add_subdirectory(csLie)
# add_subdirectory(...)
