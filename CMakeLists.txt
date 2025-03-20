cmake_minimum_required(VERSION 3.31)
project(csLieLib VERSION 1.0.0 LANGUAGES CXX)

# add the csCMake system
include(csCMake)

cs_initialize_project(csLieLib)

# prepare custom project asserts, note: use standard find_package
find_package(csCore REQUIRED COMPONENTS assert CONFIG)

list(APPEND CMAKE_MODULE_PATH "${csCore_BASE_DIR}/cmake/csCore/assert/")
include(csAssert)
csAssertGenerateProject(csLie csLie) #generate target csLie::csLieAssert
list(POP_BACK CMAKE_MODULE_PATH)


# add libs, executables, etc
add_subdirectory(csLie)
add_subdirectory(tests)
# add_subdirectory(...)

# remember to add this at the end of the main cmake file
cs_finalize_project(csLie)
