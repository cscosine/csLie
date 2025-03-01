project(csLie LANGUAGES CXX)

cs_add_library(csLie STATIC)
cs_add_library(csLie::csLie ALIAS csLie)

cs_find_package(csLie PUBLIC Eigen3 REQUIRED CONFIG)
cs_find_package(csLie PUBLIC csCore COMPONENTS csCore REQUIRED CONFIG)


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
		#Adjoint.hpp
		#Constants.hpp
		#SE2Alg.hpp
		#SE2Exp.hpp
		#SE2.hpp
		#SE2Params.hpp
		#SE3Alg.hpp
		#SE3Exp.hpp
		#SE3.hpp
		#SE3Log.hpp
		#SE3Params.hpp
		#SkewSymm.hpp
		#SO3Alg.hpp
		#SO3Exp.hpp
		#SO3.hpp
		#SO3Log.hpp
		#SO3Params.hpp
		#Utils.hpp
    #	Tagged.hpp
    S2.hpp
    #S2Alg.hpp
    #S2Exp.hpp
	PRIVATE
		#Adjoint.cpp
		#Constants.cpp
		#SE2Exp.cpp
		#SE2Params.cpp
		#SE3Exp.cpp
		#SE3Log.cpp
		#SE3Params.cpp
		#SkewSymm.cpp
		#SO3Exp.cpp
		#SO3Log.cpp
		#SO3Params.cpp
		#Utils.cpp
    S2.cpp
    #S2Exp.cpp
    )

cs_target_link_libraries(csLie
	PUBLIC
		Eigen3::Eigen
    csCore::csCore
		csLieLib::csLieLibBaseInclude
)

cs_install_target(csLie csLieExport DEFAULT)

cs_finalize_target(csLie)