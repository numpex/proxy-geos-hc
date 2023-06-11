cmake_minimum_required( VERSION 3.22.1 )

project( sem )

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------

set (CMAKE_VERBOSE_MAKEFILE ON )

#-----------------------------------------------------
# Setup the flags for users to select programming models
#-----------------------------------------------------

# use OMP to parallelise loops
set (SEM_USE_OMP OFF CACHE BOOL "")

# use RAJA to expose parallelism loop
set (SEM_USE_RAJA OFF CACHE BOOL "")

# use KOKOS to expose parallelism loop
set (SEM_USE_KOKKOS ON CACHE BOOL "")

#-----------------------------------------------------
# Setup the flags for users to select data containers
#-----------------------------------------------------

# use std::vector for vectores and arrays
set (SEM_USE_VECTOR ON CACHE BOOL "")

# use LvArray for vectors and arrays
set (SEM_USE_LVARRAY OFF CACHE BOOL "")

#-----------------------------------------------------
# Setup the flags for profiling using caliper
#-----------------------------------------------------
set (SEM_USE_CALIPER ON CACHE BOOL "")

# setup the flags for python interface
set (SEM_USE_PYTHON OFF CACHE BOOL "")

add_subdirectory(src)


