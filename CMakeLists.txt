cmake_minimum_required(VERSION 3.31)
project(csLieLib VERSION 1.0.0 LANGUAGES CXX)

# add the csCMake system
include(csCMake)

cs_initialize_project(csLieLib)

# add libs, executables, etc
add_subdirectory(csLie)
# add_subdirectory(...)

# remember to add this at the end of the main cmake file
cs_finalize_project(csLieLib)
