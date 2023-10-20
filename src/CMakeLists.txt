#-------------------------------------------------------------------
# The Solver Implementation is with Programming Models:
#
#  Users can choose from "OMP", "RAJA", and "KOKKOS" as following:
#
#  (0) The default is using sequential mode
#  (1) If SEM_USE_OMP is ON, it is using OMP to parallelise loops.
#  (2) If SEM_USE_RAJA is ON, it is using RAJA.
#  (3) If SEM_USE_KOKKOS is ON, it is using KOKKOS.
#
#  And the users can chose data containers as following:
#
#  (a) If SEM_USE_LVARRAY is ON, LVARRAY is used.
#  (b) If SEM_USE_VECTOR is ON, std::vector is used.
#
#  env variables to be defined before build
#  the variables above are defined for Pangea3
#  SEM_TPL must be first compiled and installed 
#
#  example:
#  export MODULEPATH=/data_local/appli_local/MTS/GEOSX/modulefiles/:$MODULEPATH
#  module load cmake/3.21.4 gcc/8.4.1 cuda/11.0.3 ompi/4.1.2 openblas/0.3.18 python4geosx/3.8.5-gcc-8.4.1-ompi-4.1.2
#  export SEM_TPL_ROOT_DIR=/users/j0024947/appli_src/geosx/codes/proxyApp/tpl4ProxyApp/installTPLNOCUDA
#  export SEM_TPL_ROOT_DIR=/users/j0024947/appli_src/geosx/codes/proxyApp/tpl4ProxyApp/installTPL
#  export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
#  export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
#  export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
#  export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/

#-------------------------------------------------------------------

# the solver is using sequential mode by default
set (SOLVER "SEQUENTIAL" CACHE PATH "" FORCE)

# CC options
set(CMAKE_C_COMPILER gcc-11 CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

# C++ options
set(CMAKE_CXX_COMPILER g++ CACHE PATH "")
set (CMAKE_CXX_FLAGS "-O3 -std=c++17  " CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG " CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")


# allow to install only one target
set (CMAKE_SKIP_INSTALL_ALL_DEPENDENCY "True" CACHE STRING "" FORCE)

################################
# BLT 
################################

# Default to C++17 if not set so GTest/GMock can build
if (NOT BLT_CXX_STD)
    set(BLT_CXX_STD "c++17" CACHE STRING "")
endif()

set( ENABLE_GTEST OFF CACHE STRING "" FORCE)
set( ENABLE_ALL_WARNINGS OFF CACHE STRING "" FORCE)
set( ENABLE_TESTS OFF CACHE STRING "" FORCE)
set( ENABLE_EXAMPLES OFF )
set( ENABLE_DOCS OFF )
set( ENABLE_BENCHMARKS OFF )
set( DISABLE_UNIT_TESTS ON)

if (SEM_USE_VECTOR)
    set (DATA_CONTAINER "VECTOR" CACHE STRING "" FORCE)
endif()

if(SEM_USE_OMP)
    set( ENABLE_OPENMP ON CACHE STRING "" FORCE)
    set (SOLVER "OMP" CACHE PATH "" FORCE)
    message(STATUS "BUILDING SOLVER including OMP ")
    set (SEM_USE_KOKKOS OFF CACHE BOOL "")
    set (SEM_USE_RAJA OFF CACHE BOOL "")
    set (ENABLE_CUDA OFF CACHE STRING "" FORCE)
    if(SEM_USE_LVARRAY)
      set( DATA_CONTAINER "LVARRAY" CACHE STRING "" FORCE)
      message(FATAL_ERROR "SEM_USE_OMP can  be used with SEM_USE_VECTOR=ON only")
    endif()
    set (DATA_CONTAINER "VECTOR" CACHE STRING "" FORCE)
endif()

# at this point User have to know if SEM_TPL is built with CUDA support.
if(SEM_USE_CUDA)
  set( ENABLE_CUDA ON CACHE STRING "" FORCE)
endif()
# Cuda options 
if (ENABLE_CUDA)
    set(CUDA_ROOT "/usr" CACHE PATH "")
    set(ENABLE_CUDA_NVTOOLSEXT ON CACHE BOOL "")
    set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_ROOT} CACHE PATH "")
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "")
    set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
    set(CMAKE_CUDA_ARCHITECTURES 50 CACHE STRING "")
    set(CUDA_ARCH sm_50 CACHE STRING "")
    set(CMAKE_CUDA_STANDARD 17 CACHE STRING "")
    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
    # --expt-extended-lambda --expt-relaxed-constexpr: important for compilation  are not set by default
    set(CMAKE_CUDA_FLAGS "-O3 -restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Xcompiler -std=c++17" CACHE STRING "")
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --relocatable-device-code=true" CACHE STRING "" FORCE)
    message(STATUS "CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS}")
    set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")
endif()

include( ../blt/SetupBLT.cmake )


################################
# RAJA
################################
if ( SEM_USE_RAJA)
   # disable OMP and KOKKOS 
   set (SEM_USE_OMP OFF CACHE BOOL "")
   set (SEM_USE_KOKKOS OFF CACHE BOOL "")
   # set solver
   message( " " )
   message( STATUS "-------------------RAJA SETTINGS---------------------------------------" )
   message( STATUS "-----------------------------------------------------------------------" )
   set (SOLVER "Raja" CACHE PATH "" FORCE)
   message(STATUS "BUILDING SOLVER including RAJA ")

   set (RAJA_DIR "$ENV{RAJA_DIR}")
   if(DEFINED RAJA_DIR)
       message(STATUS "SEM USE RAJA in RAJA_DIR = ${RAJA_DIR}")
       find_package(RAJA REQUIRED
                    PATHS ${RAJA_DIR}
                    NO_DEFAULT_PATH)

	    message(STATUS "SEM USE RAJA:  RAJA VERSION=${RAJA_VERSION}")

       get_target_property(RAJA_INCLUDE_DIRS RAJA INTERFACE_INCLUDE_DIRECTORIES)
       set_target_properties(RAJA PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${RAJA_INCLUDE_DIRS}")
       set(ENABLE_RAJA ON CACHE BOOL "")
       set(thirdPartyLibs ${thirdPartyLibs} RAJA )
   else()
       message(FATAL_ERROR "GEOSX requires RAJA, set RAJA_DIR to the RAJA installation directory.")
   endif()

   ################################
   # Umpire
   ################################
   set (UMPIRE_DIR "$ENV{UMPIRE_DIR}")
   if(DEFINED UMPIRE_DIR)
       message(STATUS "SEM USE RAJA: RAJA must be compiled with UMPIRE in UMPIRE_DIR = ${UMPIRE_DIR}")

       find_package(umpire REQUIRED
                    PATHS ${UMPIRE_DIR}
                    NO_DEFAULT_PATH)

	    message(STATUS "SEM USE RAJA: RAJA is compiled with UMPIRE VERSION=${umpire_VERSION}")

       set(ENABLE_UMPIRE ON CACHE BOOL "")
       set(thirdPartyLibs ${thirdPartyLibs} umpire)
   else()
       message(FATAL_ERROR "GEOSX requires Umpire, set UMPIRE_DIR to the Umpire installation directory.")
   endif()

   ################################
   # CHAI
   ################################
   set (CHAI_DIR "$ENV{CHAI_DIR}")
   if(DEFINED CHAI_DIR)
       message(STATUS "SEM USE RAJA in CHAI_DIR = ${CHAI_DIR}")

       find_package(chai REQUIRED
                    PATHS ${CHAI_DIR}
                    NO_DEFAULT_PATH)

	    message(STATUS "SEM USE RAJA: RAJA is compiled with chai_VERSION=${chai_VERSION}")

       get_target_property(CHAI_INCLUDE_DIRS chai INTERFACE_INCLUDE_DIRECTORIES)
       set_target_properties(chai PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CHAI_INCLUDE_DIRS}")

       set(ENABLE_CHAI ON CACHE BOOL "")
       set(thirdPartyLibs ${thirdPartyLibs} chai)
   else()
     message(FATAL_ERROR "GEOSX requires CHAI, set CHAI_DIR to the CHAI installation directory.")
   endif()

   if (SEM_USE_LVARRAY)
      set( DATA_CONTAINER "LVARRAY" CACHE STRING "" FORCE)
      set(ENABLE_LVARRAY ON CACHE BOOL "" FORCE)
      add_subdirectory( LvArray )
   endif()
   if (SEM_USE_VECTOR)
      set( DATA_CONTAINER "VECTOR" CACHE STRING "" FORCE)
   endif()
   message( " " )
   message( STATUS "-------------------END RAJA SETTINGS-----------------------------------" )
   message( STATUS "-----------------------------------------------------------------------" )
endif()


################################
# KOKKOS
################################
# the solver is using KOKKOS if SEM_USE_KOKKOS is ON
if (NOT SEM_USE_KOKKOS)
  message ("-- SEM_USE_KOKKOS flag is set to be false and KOKKOS is not enabled")
elseif (EXISTS $ENV{KOKKOS_DIR}/KokkosConfig.cmake)
  message ("-- SEM_USE_KOKKOS flag is set to be true and KOKOS is now enabled")
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

  if (SEM_USE_KOKKOSVECTOR)
     set (DATA_CONTAINER "KOKKOSVECTOR" CACHE PATH "" FORCE)
  endif()
  set (kokkosList "-L${Kokkos_ROOT}/../.. -lkokkoscontainers -lkokkoscore -lkokkossimd")
  set (ENABLE_KOKKOS ON CACHE BOOL "")


elseif (NOT DEFINED $ENV{KOKKOS_DIR})
  set (SEM_USE_KOKKOS OFF CACHE BOOL "")
  message("-- KOKKOS_DIR is not defined  not enabled")
else ()
  set (SEM_USE_KOKKOS OFF CACHE BOOL "")
  message("-- KOKKOS_DIR is not the correct path and is not enabled")
endif()

add_subdirectory( sem )

#add_subdirectory( tests )
