cmake_minimum_required( VERSION 3.22.1 )

project( proxyapp )

#-----------------------------------------------------
# Setup Flags for whether considering an external libs 
# for the submodule
#-----------------------------------------------------

set(USE_LVARRAY_EXT OFF CACHE BOOL "Specify if we are using an externally installed library")
set(LVARRAY_DIR "" CACHE PATH "")

#-----------------------------------------------------
# Setup Flags for Building and Installation
#-----------------------------------------------------
set (GUIX_INSTALLED_TPL OFF CACHE BOOL "Specify if we are using some TPLs installed with GUIX")
set (CMAKE_VERBOSE_MAKEFILE ON )

#-----------------------------------------------------
# Setup the flags for users to select programming models
#-----------------------------------------------------

# use OMP to parallelise loops
option (USE_OMP "Use OpenMP to parallelise loops" OFF)

# use RAJA to expose parallelism loop
option (USE_RAJA "Use RAJA to parallelise loops" OFF)

# use KOKOS to expose parallelism loop
option (USE_KOKKOS "Use KOKKOS to parallelise loops" OFF)

#---------------------------------------------------
# Tests directory to perform differen tests before
# integration into proxy apps
#---------------------------------------------------
option (ENABLE_TESTS "add Tests directory" ON)

#--------------------------------------------------------
# Option to save snapshots of wavefields (for debugging)
#--------------------------------------------------------
option (SEM_SAVE_SNAPSHOTS "save snapshots for SEM-proxy" OFF)
option (FD_SAVE_SNAPSHOTS "save snapshots for FD-proxy" OFF)

#------------------------------------------------------
# Option to printout memory allocation info of arrays
#------------------------------------------------------
option (PRINT_ALLOC_INFO "printout memory allocation info" OFF)

#---------------------------------------------------
# Option to test SEM-proxy 2D case (default is 3D)
#---------------------------------------------------
option (SEM2D "SEM-proxy 2D or 3D case" OFF)

#---------------------------------------------------------
# Option to use coloring scheme to avoid atomic operations
#---------------------------------------------------------
option (SEM_MESHCOLOR "use coloring scheme" OFF)

#---------------------------------------------------------
# Option to use hierarchical parallelism in Kokkos
#---------------------------------------------------------
option (USE_KOKKOS_TEAMS "use hierarchical parallelism in Kokkos" OFF)
enable_testing()
add_subdirectory(src)


