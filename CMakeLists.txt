cmake_minimum_required( VERSION 3.8 )

project( sem )

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "-O3 -std=c++14 " CACHE STRING "")

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------

# setup the flags for users to select programming models
set (SEM_USE_OMP ON CACHE BOOL "")
set (SEM_USE_RAJA ON CACHE BOOL "")

# setup the flags for profiling using caliper
set (SEM_USE_CALIPER ON CACHE BOOL "")

add_subdirectory(src)


