cmake_minimum_required(VERSION 3.31)
project(csLieLib VERSION 1.0.0 LANGUAGES CXX)

# add the csCMake system
include(csCMake)

cs_initialize_project(csLieLib
    # add the default configurations supported by this project
    DEFAULT_CONFIG Debug
    DEFAULT_CONFIG Release
    # and the custom configurations supported by this project, 
    # specifiying on which base configuration they are based as second value, 
    # e.g. CUSTOM_CONFIG Paranoid Debug
    CUSTOM_CONFIG Paranoid BASED_ON "Debug,Release," POSTFIX "dp"
)

# add libs, executables, etc
add_subdirectory(csLie)
# add_subdirectory(...)

# remember to add this at the end of the main cmake file
cs_finalize_project(csLie)
