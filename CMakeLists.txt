cmake_minimum_required( VERSION 3.22.1 )

project( sem )

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------

set (CMAKE_VERBOSE_MAKEFILE ON )

#-----------------------------------------------------
# architecture 
#-----------------------------------------------------
option (X86 "X86 architecture  default" ON)
option (P9 "Power9 architecture  " OFF)
#-----------------------------------------------------
# Setup the flags for users to select programming models
#-----------------------------------------------------

# use OMP to parallelise loops
option (_USE_OMP "Use OpenMP to parallelise loops" OFF)

# use RAJA to expose parallelism loop
option (_USE_RAJA "Use RAJA to parallelise loops" OFF)

# use KOKOS to expose parallelism loop
option (_USE_KOKKOS "Use KOKKOS to parallelise loops" OFF)

#-----------------------------------------------------
# Setup the flags for users to select data containers
#-----------------------------------------------------

# use std::vector for vectors and arrays
option (_USE_VECTOR "Use standard vectors as data containers" OFF)

# use LvArray for vectors and arrays
option (_USE_LVARRAY "Use LVARRAY as data containers" OFF)

#-----------------------------------------------------
# Setup the flags for profiling using caliper
#-----------------------------------------------------
option (_USE_CALIPER "Use CALIPER for profiling" OFF)

# setup the flags for python interface (not working yet)
#option (_USE_PYTHON "bind python" OFF)


add_subdirectory(src)


