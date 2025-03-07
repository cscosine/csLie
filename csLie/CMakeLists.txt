project(csLie LANGUAGES CXX)

cs_add_library(csLie INTERFACE)
cs_add_library(csLie::csLie ALIAS csLie)

cs_find_package(csLie PUBLIC Eigen3 REQUIRED CONFIG)
cs_find_package(csLie PUBLIC csCore COMPONENTS csCore REQUIRED CONFIG)
cs_find_package(csLie PRIVATE  libassert REQUIRED CONFIG)

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