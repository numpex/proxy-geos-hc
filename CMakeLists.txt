cmake_minimum_required( VERSION 3.8 )

project( sem )

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------
<<<<<<< HEAD
set( CMAKE_VERBOSE_MAKEFILE on )
=======
set (CMAKE_VERBOSE_MAKEFILE on )

>>>>>>> feature/henri/lvarray
# setup the flags for users to select programming models
# use OMP to parallelise loops
set (SEM_USE_OMP OFF CACHE BOOL "")
<<<<<<< HEAD
set (SEM_USE_RAJA ON CACHE BOOL "")
=======

# use RAJA to expose parallelism loop
set (SEM_USE_RAJA ON CACHE BOOL "")

# use std::vector for vectores and arrays
set (SEM_USE_VECTOR OFF CACHE BOOL "")

#use LvArray for vectors and arrays
set (SEM_USE_LVARRAY ON CACHE BOOL "")
>>>>>>> feature/henri/lvarray

# setup the flags for profiling using caliper
set (SEM_USE_CALIPER ON CACHE BOOL "")

add_subdirectory(src)


