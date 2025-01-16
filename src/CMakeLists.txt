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

if( ENABLE_CUDA AND NOT USE_KOKKOS)
  set(extraComponentsLinkList ${extraComponentsLinkList} cuda )
  message ("-- ENABLE_CUDA flag is set: extraComponentsLinkList= ${extraComponentsLinkList} ")
endif()

include( ../blt/SetupBLT.cmake )


################################
# SEM USE SEQUENTIAL MODEL
################################
# the solver is sequential +Vector on CPU

if(NOT (USE_OMP OR USE_KOKKOS OR USE_RAJA))
	message(STATUS "BUILDING SEQUENTIAL SOLVER including Vector on CPU")
	if(ENABLE_OPENMP OR ENABLE_CUDA OR ENABLE_MPI OR ENABLE_HIP)
		message(FATAL_ERROR "The sequential mode is not compatible with any GPU-accelerated options or parallel model")
	endif()
endif()

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
  set (USE_VECTOR ON CACHE BOOL "" FORCE)

  # We could avoid forcing the value of ENABLE_OPENMP to consider the purely sequential code ?
  #set( ENABLE_OPENMP ON CACHE STRING "" FORCE)

  if(ENABLE_CUDA OR ENABLE_HIP OR NOT ENABLE_OPENMP)
	  message("-- USE_OMP is not compatible with GPU acceleration are going to be turned OFF: ENABLE_CUDA = " ${ENABLE_CUDA} " ENABLE_HIP = ")
	  set (ENABLE_CUDA OFF CACHE STRING "" FORCE)
	  set (ENABLE_HIP OFF CACHE STRING "" FORCE)
	  set( ENABLE_OPENMP ON CACHE STRING "" FORCE)
	  
	  
	  # Include the core config file after all the option are set
	  message("-- Including the config file to update the compilers args according to the forced options")
	  ## The directory with the config. files
	  set(CONFIG_DIR ${CMAKE_SOURCE_DIR}/../$ENV{proxy_tpl_repo}/configs)
	  if(EXISTS ${CONFIG_DIR}/config_core.cmake)
		  include(${CONFIG_DIR}/config_core.cmake)
	  else()
		  message(FATAL_ERROR "The config_core.cmake is not found in the the config files directory: " ${CONFIG_DIR})
	  endif()
  endif()
  #set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
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
   if(NOT (${ENABLE_OPENMP} OR ${ENABLE_CUDA}))
           message(STATUS "ENABLE_OPENMP is set to " ${ENABLE_OPENMP} " ENABLE_CUDA " ${ENABLE_CUDA})
	   message(FATAL_ERROR "The TPLs have been built without enabling one of the OPENMP or CUDA feature required for USE_RAJA option: please rebuild the RAJA library.")
   endif()
   if(DEFINED RAJA_DIR)
       message(STATUS "USE RAJA in RAJA_DIR = ${RAJA_DIR}")
       find_package(RAJA REQUIRED
                    PATHS ${RAJA_DIR}
                    NO_DEFAULT_PATH)
       get_target_property(RAJA_INCLUDE_DIRS RAJA INTERFACE_INCLUDE_DIRECTORIES)
       set_target_properties(RAJA PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${RAJA_INCLUDE_DIRS}")
       message("-- RAJA_INCLUDE_DIRS: ${RAJA_INCLUDE_DIRS}")
       # To get the link libraries of RAJA: blt::cuda, blt::openmp, RAJA::cub; camp; dl
       #get_target_property(RAJA_INTERFACE_LINK_LIBRARIES RAJA INTERFACE_LINK_LIBRARIES)
       #message(STATUS "RAJA_INTERFACE_LINK_LIBRARIES " ${RAJA_INTERFACE_LINK_LIBRARIES})
       #set (dependencyList ${dependencyList} ${RAJA_INTERFACE_LINK_LIBRARIES} )   
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

   #set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
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
elseif (DEFINED KOKKOS_DIR)
  set (USE_LVARRAY OFF CACHE BOOL "" FORCE)
  set (USE_VECTOR OFF CACHE BOOL "" FORCE)
  message ("-- USE_KOKKOS flag is set to be true and KOKOS is now enabled")
  message("-- Found the KOKKOS library in ${KOKKOS_DIR}")
  set (SOLVER "Kokkos" CACHE PATH "" FORCE)

  # KOKOS uses OMP
  #find_package(OpenMP) is not required regarding the previous set up of BLT
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  message(STATUS "USE KOKKOS in KOKKOS_DIR = ${KOKKOS_DIR}")
  
  # setup KOKOS
  find_package(Kokkos REQUIRED
                    PATHS ${KOKKOS_DIR}
                    NO_DEFAULT_PATH)
  get_target_property(KOKKOS_INCLUDE_DIRS Kokkos::kokkos INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "KOKKOS_INCLUDE_DIRS: ${KOKKOS_INCLUDE_DIRS}")     


  get_target_property(KOKKOS_LINK_LIBS Kokkos::kokkos INTERFACE_LINK_LIBRARIES)
  message("-- KOKKOS_LINK_LIBS: ${KOKKOS_LINK_LIBS}")

  #get_target_property(KOKKOSCORE_LINK_OPTION Kokkos::kokkoscore INTERFACE_LINK_OPTIONS)
  #message(STATUS "KOKKOSCORE_LINK_OPTION: ${KOKKOSCORE_LINK_OPTION}")

  set (extraComponentsLinkList ${extraComponentsLinkList} Kokkos::kokkos)
  #set(thirdPartyLibs ${thirdPartyLibs} Kokkos::kokkos )

  ## -- The idea below is to get the self-contained dependencies of Kokkos (Kokkos::openmp, Kokkos::cuda etc..)
  #get_target_property(KOKKOSCORE_LINK_LIBS Kokkos::kokkoscore INTERFACE_LINK_LIBRARIES)
  #message("-- KOKKOSCORE_LINK_LIBS" ${KOKKOSCORE_LINK_LIBS})
  #set (dependencyList ${dependencyList} ${KOKKOSCORE_LINK_LIBS})



  #set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
  message (STATUS "KOKKOS DEPENDENCY LIST= ${dependencyList}")
  if( ENABLE_HIP )
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
    set (dependencyList ${dependencyList} Kokkos::kokkos)
    message (STATUS "KOKKOS DEPENDENCY LIST= ${dependencyList}")
    include_directories("/opt/rocm-6.0.0/include")
  endif()

else ()
  set (USE_KOKKOS OFF CACHE BOOL "")
  message("-- KOKKOS_DIR is not the correct path and is not enabled / KOKKOS_DIR = " ${KOKKOS_DIR})
endif()

message(STATUS "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

configure_file( utils/commonConfig.hpp.in
                ${CMAKE_BINARY_DIR}/src/commonConfig.hpp)

add_subdirectory( sem )
add_subdirectory( fd )
if ( addTests)
  add_subdirectory( tests )
endif()
