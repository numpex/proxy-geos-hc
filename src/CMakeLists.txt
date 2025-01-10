#-------------------------------------------------------------------
# The Solver Implementation is with Programming Models:
#
#  Users can choose from "OMP", "RAJA", and "KOKKOS" as following:
#
#  (0) The default is using sequential mode
#  (1) If USE_OMP is ON, it is using OpenMP+Vector on CPU.
#  (2) If USE_RAJA is ON, it is using RAJA+LVARRAY+CUDA.
#  (3) If USE_KOKKOS is ON, it is using KOKKOS+CUDA.
#
#  And the users can chose data containers as following:
#
#  env variables to be defined before build
#  the variables above are defined for Pangea3
#  _TPL must be first compiled and installed
#
#  example:
#  export MODULEPATH=/data_local/appli_local/MTS/GEOSX/modulefiles/:$MODULEPATH
#  module load cmake/3.21.4 gcc/8.4.1 cuda/11.0.3 ompi/4.1.2 openblas/0.3.18 python4geosx/3.8.5-gcc-8.4.1-ompi-4.1.2
#  export _TPL_ROOT_DIR=/users/j0024947/appli_src/geosx/codes/proxyApp/tpl4ProxyApp/installTPL
#  export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
#  export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
#  export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
#  export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
#
#-------------------------------------------------------------------

# the solver is using sequential mode by default
set (SOLVER "SEQUENTIAL" CACHE PATH "" FORCE)
set (USE_VECTOR ON CACHE BOOL "" FORCE)

# allow to install only one target
set (CMAKE_SKIP_INSTALL_ALL_DEPENDENCY "True" CACHE STRING "" FORCE)

################################
# BLT
################################

# Default to C++17 if not set so GTest/GMock can build
if (NOT BLT_CXX_STD)
    set(BLT_CXX_STD "c++17" CACHE STRING "")
endif()

# the following options are setup for BLT
set( ENABLE_GTEST OFF CACHE STRING "" FORCE)
set( ENABLE_ALL_WARNINGS OFF CACHE STRING "" FORCE)
set( ENABLE_TESTS OFF CACHE STRING "" FORCE)
set( ENABLE_EXAMPLES OFF )
set( ENABLE_DOCS OFF )
set( ENABLE_BENCHMARKS OFF )
set( DISABLE_UNIT_TESTS ON)
# force ENABLE_OPENMP to ON
set( ENABLE_OPENMP ON CACHE STRING "" FORCE)

include( ../blt/SetupBLT.cmake )

################################
# SEM USE OMP
################################
# the solver is using OpenMP+Vector on CPU
if(USE_OMP)
  set (SOLVER "OMP" CACHE PATH "" FORCE)
  message(STATUS "BUILDING SOLVER including OpenMP+Vector on CPU")
  set (USE_KOKKOS OFF CACHE BOOL "")
  set (USE_RAJA OFF CACHE BOOL "")
  set (USE_LVARRAY OFF CACHE BOOL "" FORCE)
  set (ENABLE_CUDA OFF CACHE STRING "" FORCE)
  set (USE_VECTOR ON CACHE BOOL "" FORCE)
  set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
  #message ("-- USE_OMP flag is set: USE_VECTOR= ${USE_VECTOR} ")
  message ("-- USE_OMP flag is set: extraComponentsLinkList= ${extraComponentsLinkList} ")
else()
    message ("-- USE_OMP flag is set to be false and OMP is not enabled")
endif()


################################
# RAJA
################################
# the solver is using RAJA+LVARRAY+CUDA on GPU
if ( USE_RAJA)
   # disable OMP and KOKKOS
   set (USE_LVARRAY ON CACHE BOOL "" FORCE)
   set (USE_VECTOR OFF CACHE BOOL "" FORCE)
   set (USE_OMP OFF CACHE BOOL "")
   set (USE_KOKKOS OFF CACHE BOOL "")
   # set solver
   message( " " )
   message( STATUS "-------------------RAJA SETTINGS---------------------------------------" )
   set (SOLVER "Raja" CACHE PATH "" FORCE)
   message(STATUS "BUILDING SOLVER including RAJA ")

   #set (RAJA_DIR "$ENV{RAJA_DIR}")
   if(DEFINED RAJA_DIR)
       message(STATUS "USE RAJA in RAJA_DIR = ${RAJA_DIR}")
       find_package(RAJA REQUIRED
                    PATHS ${RAJA_DIR}
                    NO_DEFAULT_PATH)
       get_target_property(RAJA_INCLUDE_DIRS RAJA INTERFACE_INCLUDE_DIRECTORIES)
       set_target_properties(RAJA PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${RAJA_INCLUDE_DIRS}")
       message(STATUS "RAJA_INCLUDE_DIRS: ${RAJA_INCLUDE_DIRS}")
       set(ENABLE_RAJA ON CACHE BOOL "")
       set(thirdPartyLibs ${thirdPartyLibs} RAJA )
   else()
       message(FATAL_ERROR "ProxyApp requires RAJA, set RAJA_DIR to the RAJA installation directory.")
   endif()

   if( ENABLE_HIP )
     set(RAJA_ENABLE_HIP ON CACHE BOOL "" FORCE)
     set(RAJA_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
     include_directories("/opt/rocm-6.0.0/include")
   endif()

   ################################
   # Umpire
   ################################
   #set (UMPIRE_DIR "$ENV{UMPIRE_DIR}")
   if(DEFINED UMPIRE_DIR)
       message(STATUS " USE RAJA: RAJA must be compiled with UMPIRE in UMPIRE_DIR = ${UMPIRE_DIR}")
       find_package(umpire REQUIRED
                    PATHS ${UMPIRE_DIR}
                    NO_DEFAULT_PATH)
       set(ENABLE_UMPIRE ON CACHE BOOL "")
       set(thirdPartyLibs ${thirdPartyLibs} umpire)
   else()
       message(FATAL_ERROR "ProxyApp requires Umpire, set UMPIRE_DIR to the Umpire installation directory.")
   endif()

   ################################
   # CHAI
   ################################
   #set (CHAI_DIR "$ENV{CHAI_DIR}")
   if(DEFINED CHAI_DIR)
       message(STATUS " USE RAJA in CHAI_DIR = ${CHAI_DIR}")
       find_package(chai REQUIRED
                    PATHS ${CHAI_DIR}
                    NO_DEFAULT_PATH)
       get_target_property(CHAI_INCLUDE_DIRS chai INTERFACE_INCLUDE_DIRECTORIES)
       set_target_properties(chai PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CHAI_INCLUDE_DIRS}")
       set(ENABLE_CHAI ON CACHE BOOL "")
       set(thirdPartyLibs ${thirdPartyLibs} chai)
   else()
     message(FATAL_ERROR "ProxyApp requires CHAI, set CHAI_DIR to the CHAI installation directory.")
   endif()

   set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
   set (dependencyList lvarray RAJA chai)
   
   ################################
   # CALIPER
   ################################
if(ENABLE_CALIPER)
    if(NOT EXISTS ${CALIPER_DIR})
        message(FATAL_ERROR "CALIPER_DIR must be defined and point to a valid directory when using caliper.")
    endif()

    message(STATUS "Using caliper from ${CALIPER_DIR}")

    find_package(caliper REQUIRED
                 PATHS ${CALIPER_DIR})
	 message(STATUS "Caliper found in CMAkeList")
    blt_register_library(NAME caliper
                         INCLUDES ${caliper_INCLUDE_PATH}
                         LIBRARIES caliper
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(thirdPartyLibs ${thirdPartyLibs} caliper)
else()
    message(STATUS "Not using caliper.")
endif()

   
   add_subdirectory( LvArray )
   include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LvArray/src)

   if( ENABLE_HIP )
     set (dependencyList ${dependencyList} blt::hip blt::hip_runtime)
   endif()

   message( STATUS "-------------------END RAJA SETTINGS-----------------------------------" )
   message( " " )
else()
    message ("-- USE_RAJA flag is set to be false and RAJA is not enabled")
endif()


################################
# KOKKOS
################################
# the solver is using KOKKOS+CUDA if USE_KOKKOS is ON
if (NOT USE_KOKKOS)
  message ("-- USE_KOKKOS flag is set to be false and KOKKOS is not enabled")
elseif (EXISTS $ENV{KOKKOS_DIR}/KokkosConfig.cmake)
  set (USE_LVARRAY OFF CACHE BOOL "" FORCE)
  set (USE_VECTOR OFF CACHE BOOL "" FORCE)
  message ("-- USE_KOKKOS flag is set to be true and KOKOS is now enabled")
  message("-- Found the KOKKOS library in $ENV{KOKKOS_DIR}")
  set (SOLVER "Kokkos" CACHE PATH "" FORCE)

  # KOKOS uses OMP
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

  # setup KOKOS
  set(Kokkos_ROOT "$ENV{KOKKOS_DIR}")
  find_package(Kokkos REQUIRED)
  set (Kokkos_INCLUDE_DIRS "${Kokkos_ROOT}/../../../include")
  include_directories(${Kokkos_INCLUDE_DIRS})

  set (kokkosList "-L${Kokkos_ROOT}/../.. -lkokkoscontainers -lkokkoscore -lkokkossimd")
  set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
  set (dependencyList ${dependencyList} ${kokkosList})
  message (STATUS "KOKKOS DEPENDENCY LIST= ${dependencyList}")
  if( ENABLE_HIP )
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
    set (dependencyList ${dependencyList} Kokkos::kokkos)
    message (STATUS "KOKKOS DEPENDENCY LIST= ${dependencyList}")
    include_directories("/opt/rocm-6.0.0/include")
  elseif ( ENABLE_CUDA )
    set (dependencyList "${dependencyList} -lcuda")
  endif()

elseif (NOT DEFINED $ENV{KOKKOS_DIR})
  set (USE_KOKKOS OFF CACHE BOOL "")
  message("-- KOKKOS_DIR is not defined  not enabled")
else ()
  set (USE_KOKKOS OFF CACHE BOOL "")
  message("-- KOKKOS_DIR is not the correct path and is not enabled")
endif()

message(STATUS "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

configure_file( utils/commonConfig.hpp.in
                ${CMAKE_BINARY_DIR}/src/commonConfig.hpp)

add_subdirectory( sem )
add_subdirectory( fd )
if ( addTests)
  add_subdirectory( tests )
endif()
