cmake_minimum_required( VERSION 3.22.1 )

project( sem )

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------

set (CMAKE_VERBOSE_MAKEFILE ON )

#set (CMAKE_C_COMPILER "/hrtc/apps/devtools/llvm/x86_64/centos8/17.0.0/gcc12/bin/clang")
#set (CMAKE_CXX_COMPILER "/hrtc/apps/devtools/llvm/x86_64/centos8/17.0.0/gcc12/bin/clang++")

#-----------------------------------------------------
# Setup the flags for users to select programming models
#-----------------------------------------------------

# use OMP to parallelise loops
option (SEM_USE_OMP "Use OpenMP to parallelise loops" OFF)

# use RAJA to expose parallelism loop
option (SEM_USE_RAJA "Use RAJA to parallelise loops" OFF)

# use KOKOS to expose parallelism loop
option (SEM_USE_KOKKOS "Use KOKKOS to parallelise loops" OFF)

#-----------------------------------------------------
# Setup the flags for users to select data containers
#-----------------------------------------------------

# use std::vector for vectores and arrays
option (SEM_USE_VECTOR "Use standard vectors as data containers" OFF)

# use LvArray for vectors and arrays
option (SEM_USE_LVARRAY "Use LVARRAY as data containers" OFF)

#-----------------------------------------------------
# Setup the flags for profiling using caliper
#-----------------------------------------------------
set (SEM_USE_CALIPER ON CACHE BOOL "")

# setup the flags for python interface
set (SEM_USE_PYTHON OFF CACHE BOOL "")

add_subdirectory(src)


