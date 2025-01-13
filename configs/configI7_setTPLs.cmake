# Set up the tpls
set( GEOSX_TPL_ROOT_DIR ${CMAKE_SOURCE_DIR}/../$ENV{proxy_tpl_repo} CACHE PATH "")

# Include the cache file of the third-party library
include(${GEOSX_TPL_ROOT_DIR}/configs/configI7.cmake) 

set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/installTPLs CACHE PATH "")

set(CAMP_DIR ${GEOSX_TPL_DIR}/raja CACHE PATH "")
set(RAJA_DIR ${GEOSX_TPL_DIR}/raja CACHE PATH "")
set( RAJA_ENABLE_VECTORIZATION OFF CACHE BOOL "" FORCE)

#set(ENABLE_UMPIRE ON CACHE BOOL "")
set(UMPIRE_DIR ${GEOSX_TPL_DIR}/chai CACHE PATH "")

#set(ENABLE_CHAI ON CACHE BOOL "")
set(CHAI_DIR ${GEOSX_TPL_DIR}/chai CACHE PATH "")

#set(ENABLE_CALIPER ON CACHE BOOL "")
#set(ENABLE_ADIAK ON CACHE BOOL "" )
set(CALIPER_DIR ${GEOSX_TPL_DIR}/caliper CACHE PATH "")
set(adiak_DIR ${GEOSX_TPL_DIR}/adiak/lib/cmake/adiak/ CACHE PATH "")
#set(adiak_DIR ${GEOSX_TPL_DIR}/adiak/lib CACHE PATH "")

# Set Kokkos_ROOT and KOKKOS_DIR
set(KOKKOS_DIR ${GEOSX_TPL_DIR}/kokkos CACHE PATH "")
message(STATUS "-- The KOKKOS_DIR is " ${KOKKOS_DIR})
# BLT as submodule can be tricky in the sens that the CUDA_ARCH need to be passed as argument ./blt/cmake/thirdparty/SetupHIP.cmake:15:          $ENV{ROCM_DIR}
# ./blt/cmake/BLTOptions.cmake:61:set(BLT_CLANG_CUDA_ARCH "$ENV{CUDA_ARCH}" CACHE STRING "Compute architecture to use when generating CUDA code with Clang")


#set(ENABLE_ADDR2LINE ON CACHE BOOL "")

#set(CHAI_CUDA_FLAGS "-arch ${CUDA_ARCH}" CACHE STRING "" FORCE)


# GTEST options
#set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")
#set(gtest_disable_pthreads ON CACHE BOOL "")

# Documentation
#set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "" FORCE)
#set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
