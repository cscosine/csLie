project(csLie LANGUAGES CXX)

cs_add_library(csLie STATIC)
cs_add_library(csLie::csLie ALIAS csLie)

cs_find_package(csLie PUBLIC Eigen3 REQUIRED CONFIG)


# use the target_sources_wrapper, which wrap cmake native target_sources. 
# Wrapping is necessary to cover both build and installation step
# Normally headers are PUBLIC and cpp are PRIVATE, but you might have private headers too
# If the option NO_AUTO_SOURCE_GROUP is omitted the target sources and headers are 
# automatically grouped respectively in "Sources" and "Headers" folder in the IDE
# environment whereas the dependency headers (i.e. those from linked targets) are grouped
# in the "Dependency Headers" folder.
# Specify NO_AUTO_SOURCE_GROUP if you want all sources and headers grouped in the IDE
# standard folders ("Source Files" and "Header Files" in Visual Studio).
# note: you can repeat target_sources_wrapper if you want to organize your files
cs_target_sources(csLie 
	PRIVATE_INSTALL 
		Adjoint.h
		Constants.h
		SE2Alg.h
		SE2Exp.h
		SE2.h
		SE2Params.h
		SE3Alg.h
		SE3Exp.h
		SE3.h
		SE3Log.h
		SE3Params.h
		SkewSymm.h
		SO3Alg.h
		SO3Exp.h
		SO3.h
		SO3Log.h
		SO3Params.h
		Utils.h
    Tagged.h
    S2.h
    S2Alg.h
    S2Exp.h
	PRIVATE
		Adjoint.cpp
		Constants.cpp
		SE2Exp.cpp
		SE2Params.cpp
		SE3Exp.cpp
		SE3Log.cpp
		SE3Params.cpp
		SkewSymm.cpp
		SO3Exp.cpp
		SO3Log.cpp
		SO3Params.cpp
		Utils.cpp
    S2.cpp
    S2Exp.cpp
    )

cs_target_link_libraries(csLie
	PUBLIC
		Eigen3::Eigen
)

add_subdirectory(tests)


cs_install_target(csLie csLieExport DEFAULT)

cs_finalize_target(csLie)