# Set up the TPLs
set( _TPL_ROOT_DIR $ENV{proxy_tpl_dir} CACHE PATH "")
if(NOT EXISTS ${_TPL_ROOT_DIR})
	message(FATAL_ERROR "The path provided for the TPLs is not valid: _TPL_ROOT_DIR = " ${_TPL_ROOT_DIR})
endif()

# Include the cache file of the third-party libraries
set(config_tpl $ENV{config_tpl} CACHE PATH "")
if(EXISTS ${_TPL_ROOT_DIR}/configs/${config_tpl})
	include(${_TPL_ROOT_DIR}/configs/${config_tpl}) 
else()
	message(FATAL_ERROR "The config_tpl file ${config_tpl} is not found in the provided _TPL_ROOT_DIR " ${_TPL_ROOT_DIR})
endif()

# To keep track of change: Beaware that _TPL_INSTALL_DIR was originally called GEOSX_TPL_DIR 
# and used in some available configs provided in the LvArray submodule (src/LvArray/cmake/blt/host-configs/)
set(_TPL_INSTALL_DIR ${_TPL_ROOT_DIR}/$ENV{install_tpl} CACHE PATH "")

set(CAMP_DIR ${_TPL_INSTALL_DIR}/raja CACHE PATH "")
set(RAJA_DIR ${_TPL_INSTALL_DIR}/raja CACHE PATH "")
set( RAJA_ENABLE_VECTORIZATION OFF CACHE BOOL "" FORCE)

#set(ENABLE_UMPIRE ON CACHE BOOL "")
set(UMPIRE_DIR ${_TPL_INSTALL_DIR}/chai CACHE PATH "")

#set(ENABLE_CHAI ON CACHE BOOL "")
set(CHAI_DIR ${_TPL_INSTALL_DIR}/chai CACHE PATH "")

#set(ENABLE_CALIPER ON CACHE BOOL "")
set(CALIPER_DIR ${_TPL_INSTALL_DIR}/caliper CACHE PATH "")

#set(ENABLE_ADIAK ON CACHE BOOL "" )
set(adiak_DIR ${_TPL_INSTALL_DIR}/adiak/lib/cmake/adiak/ CACHE PATH "")

# Set Kokkos_ROOT and KOKKOS_DIR
set(KOKKOS_DIR ${_TPL_INSTALL_DIR}/kokkos CACHE PATH "")

#set(ENABLE_ADDR2LINE ON CACHE BOOL "")

#set(CHAI_CUDA_FLAGS "-arch ${CUDA_ARCH}" CACHE STRING "" FORCE)

#ENABLE_FIND_MPI to control the setting of the MPI related compiler flag and options. 
#Be aware that "FindMPI can break nvcc. In that case, you should set ENABLE_FIND_MPI to Off and control" (./cmake/thirdparty/SetupCUDA.cmake #87)


# GTEST options
#set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")
#set(gtest_disable_pthreads ON CACHE BOOL "")

# Documentation
#set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "" FORCE)
#set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
