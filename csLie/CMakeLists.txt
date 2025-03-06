project(csLie LANGUAGES CXX)

cs_add_library(csLie INTERFACE)
cs_add_library(csLie::csLie ALIAS csLie)

cs_find_package(csLie PUBLIC Eigen3 REQUIRED CONFIG)
cs_find_package(csLie PUBLIC csCore COMPONENTS csCore REQUIRED CONFIG)
cs_find_package(csLie PRIVATE  libassert REQUIRED CONFIG)


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
    	csCore::csCore
		csLieLib::csLieLib_csCMakeBaseInclude
    	libassert::assert
)

cs_install_target(csLie csLieExport DEFAULT)

cs_finalize_target(csLie)