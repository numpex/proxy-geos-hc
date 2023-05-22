cmake_minimum_required( VERSION 3.8 )

project( sem )

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------
set (CMAKE_VERBOSE_MAKEFILE on )

# setup the flags for users to select programming models
set (SEM_USE_OMP OFF CACHE BOOL "")
set (SEM_USE_RAJA ON CACHE BOOL "")

# setup the flags for profiling using caliper
set (SEM_USE_CALIPER ON CACHE BOOL "")

add_subdirectory(src)


