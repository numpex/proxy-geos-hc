cmake_minimum_required( VERSION 3.8 )

project( sem )

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------
set (CMAKE_VERBOSE_MAKEFILE on )

# setup the flags for users to select programming models
# use OMP to parallelise loops
set (SEM_USE_OMP ON CACHE BOOL "")

# use RAJA to expose parallelism loop
set (SEM_USE_RAJA OFF CACHE BOOL "")

# use std::vector for vectores and arrays
set (SEM_USE_VECTOR ON CACHE BOOL "")

#use LvArray for vectors and arrays
set (SEM_USE_LVARRAY OFF CACHE BOOL "")

# setup the flags for profiling using caliper
set (SEM_USE_CALIPER ON CACHE BOOL "")

# setup the flags for python interface
set (SEM_USE_PYTHON OFF CACHE BOOL "")

add_subdirectory(src)


