cmake_minimum_required( VERSION 3.22.1 )

project( proxyapp)

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------

set (CMAKE_VERBOSE_MAKEFILE ON )

#-----------------------------------------------------
# Setup the computing resources for running proxyapp
#-----------------------------------------------------
option (X86_laptop "X86 architecture  default" OFF)
option (X86_cypress "Compilation on Cypress Cluster" ON)
option (P9 "Power9 architecture  " OFF)
option (ARM "ARM aarch64 architecture: e.g. grace CPU" OFF)

#-----------------------------------------------------
# Setup the flags for users to select programming models
#-----------------------------------------------------

# use OMP to parallelise loops
option (USE_OMP "Use OpenMP to parallelise loops" OFF)

# use RAJA to expose parallelism loop
option (USE_RAJA "Use RAJA to parallelise loops" OFF)

# use KOKOS to expose parallelism loop
option (USE_KOKKOS "Use KOKKOS to parallelise loops" OFF)

#-----------------------------------------------------
# Setup the flags for profiling using caliper
#-----------------------------------------------------
option (USE_CALIPER "Use CALIPER for profiling" OFF)

#---------------------------------------------------
# Tests directory to perform differen tests before 
# integration into proxy apps
#---------------------------------------------------
option (addTests " add Tests directory" OFF)
option (addFD " add finite difference proxy directory" OFF)


add_subdirectory(src)


