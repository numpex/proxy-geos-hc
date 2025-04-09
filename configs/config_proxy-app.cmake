# Include the config file for pre-loading the cache
if(DEFINED ENV{config_proxy})
	set(config_proxy $ENV{config_proxy} CACHE PATH "")
	if(EXISTS ${config_proxy})
		include(${config_proxy})
	else()	
		message(FATAL_ERROR "The config_proxy file ${config_proxy} is not found: Please check the provided path")
	endif()
else()
	message(STATUS "---No config file specified--- = "${config_proxy} " or "$ENV{config_proxy})
endif()
#####################################
############# What is enabled for RAJA
######################################
if(USE_RAJA)
	set( RAJA_ENABLE_VECTORIZATION OFF CACHE BOOL "" FORCE)
	set(ENABLE_UMPIRE ON CACHE BOOL "" FORCE)
	set(ENABLE_CHAI ON CACHE BOOL "" FORCE)
	set(ENABLE_CALIPER OFF CACHE BOOL "" FORCE)
	set(ENABLE_ADIAK OFF CACHE BOOL "" FORCE)
	
	#Inherited from the TPLs config file   
	#set(RAJA_ENABLE_CUDA ${ENABLE_CUDA} CACHE BOOL "" FORCE)
	#set(RAJA_ENABLE_OPENMP ${ENABLE_OPENMP} CACHE BOOL "" FORCE)
endif()
message(STATUS "BUILD_FROM_TPLMIRROR " ${BUILD_FROM_TPLMIRROR})
if(BUILD_FROM_TPLMIRROR OR DEFINED ENV{install_tpl_dir})
	message(STATUS "--Setting the paths for the TPL")
	# To keep track of change: Beaware that _TPL_INSTALL_DIR was originally called GEOSX_TPL_DIR 
	# and used in some available configs provided in the LvArray submodule (src/LvArray/cmake/blt/host-configs/)
	set(_TPL_INSTALL_DIR $ENV{install_tpl_dir} CACHE PATH "")
	
	if(USE_RAJA)
		set(CAMP_DIR ${_TPL_INSTALL_DIR}/raja CACHE PATH "")
		set(RAJA_DIR ${_TPL_INSTALL_DIR}/raja CACHE PATH "")
		
		if(ENABLE_UMPIRE OR ENABLE_CHAI)
			set(UMPIRE_DIR ${_TPL_INSTALL_DIR}/chai CACHE PATH "")
		endif()
		
		if(ENABLE_CHAI)
			#set(CHAI_DIR ${_TPL_INSTALL_DIR}/chai/share/chai/cmake CACHE PATH "")
			set(CHAI_DIR ${_TPL_INSTALL_DIR}/chai CACHE PATH "")
		endif()
	
		if(ENABLE_CALIPER)
			set(CALIPER_DIR ${_TPL_INSTALL_DIR}/caliper CACHE PATH "")
		endif()
		
		#if(ENABLE_ADIAK)
			set(adiak_DIR ${_TPL_INSTALL_DIR}/adiak/lib/cmake/adiak/ CACHE PATH "")
			#endif()
		set(ENABLE_ADDR2LINE ON CACHE BOOL "")
	endif()
	
	if(USE_KOKKOS)
		# Set Kokkos_ROOT and KOKKOS_DIR
		set(KOKKOS_DIR ${_TPL_INSTALL_DIR}/kokkos CACHE PATH "")
	endif()
else()
	message(STATUS "Guix-installed TPLs: the paths are not set for find_package")
endif()

#set(CHAI_CUDA_FLAGS "-arch ${CUDA_ARCH}" CACHE STRING "" FORCE)

#ENABLE_FIND_MPI to control the setting of the MPI related compiler flag and options. 
#Be aware that "FindMPI can break nvcc. In that case, you should set ENABLE_FIND_MPI to Off and control" (./cmake/thirdparty/SetupCUDA.cmake #87)


# GTEST options
#set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")
#set(gtest_disable_pthreads ON CACHE BOOL "")

# Documentation
#set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "" FORCE)
#set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
