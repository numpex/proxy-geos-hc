#-----------------------------------------------------------------------------------------------------------------------------
# The Solver Implementation is with Programming Models:
#
#  Users can choose from "OMP", "RAJA", and "KOKKOS" as following:
#
#  (0) The default is using sequential mode
#  (1) If USE_OMP is ON, it is using OpenMP+Vector on CPU.
#  (2) If USE_RAJA is ON, it is using RAJA+LVARRAY+CUDA.
#  (3) If USE_KOKKOS is ON, it is using KOKKOS with or without CUDA and OpenMP.
#
#
#  Edit the env_var.sh file in the TPLs folder to set the Env variables to be defined before build of both TPL and the ProxyApp 
#  _TPL must be first compiled and installed from the TPLs repo
#  The exported Env variables are required to set the _TPL paths in the config file configs/config_proxy-app.cmake 
#
#-----------------------------------------------------------------------------------------------------------------------------

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

if(USE_OMP)
  if(NOT ENABLE_OPENMP)
    set(COPY_ENABLE_OPENMP ${ENABLE_OPENMP} CACHE BOOL "" FORCE)
    set(ENABLE_OPENMP ON CACHE BOOL "Forcing ENABLE_OPENMP before BLT Setup for the option USE_OMP" FORCE)
    message(STATUS "--Forcing ENABLE_OPENMP prior to BLT setup")
  endif()
endif()

if (DEFINED BLT_SOURCE_DIR)
    # Support having a shared BLT outside of the repository if given a BLT_SOURCE_DIR

    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake: " ${BLT_SOURCE_DIR})
    endif()
else()
    # Use internal BLT if no BLT_SOURCE_DIR is given
    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR
            "The BLT submodule is not present. "
            "Run the following two commands in your git repository: \n"
            "    git submodule init\n"
            "    git submodule update" )
    endif()
endif()

include(${BLT_SOURCE_DIR}/SetupBLT.cmake )

set(dependencyList "")
set(extraComponentsLinkList "")

if(DEFINED COPY_ENABLE_OPENMP)
	set(ENABLE_OPENMP ${COPY_ENABLE_OPENMP} CACHE BOOL "Set back ENABLE_OPENMP to its original OFF value" FORCE)
endif()
message(STATUS "The cuda RT library is " ${CUDA_LIBRARIES})

if( ENABLE_CUDA AND NOT USE_KOKKOS)
  set(extraComponentsLinkList ${extraComponentsLinkList} cuda )
  message ("-- ENABLE_CUDA flag is set: extraComponentsLinkList= ${extraComponentsLinkList} ")
endif()
################################
# SEM USE SEQUENTIAL MODEL
################################
# the solver is sequential +Vector on CPU

if(NOT (USE_OMP OR USE_KOKKOS OR USE_RAJA))
	message(STATUS "BUILDING SEQUENTIAL SOLVER including Vector on CPU")
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

  set (dependencyList ${dependencyList} openmp )
  message ("-- USE_OMP flag is set: dependencyList= ${dependencyList} ")
else()
  message ("-- USE_OMP flag is set to be false and OMP is not enabled")
endif()

################################
# RAJA
################################
# the solver is using RAJA+LVARRAY+CUDA on GPU
if ( USE_RAJA)
   set (USE_LVARRAY ON CACHE BOOL "" FORCE)
   set (USE_VECTOR OFF CACHE BOOL "" FORCE)
   # disable OMP and KOKKOS
   set (USE_OMP OFF CACHE BOOL "")
   set (USE_KOKKOS OFF CACHE BOOL "")
  
   # set solver
   message( " " )
   message( STATUS "-------------------RAJA SETTINGS---------------------------------------" )
   set (SOLVER "Raja" CACHE PATH "" FORCE)
   message(STATUS "BUILDING SOLVER including RAJA with the following setting: RAJA_ENABLE_CUDA = "${RAJA_ENABLE_CUDA} " and RAJA_ENABLE_OPENMP = " ${RAJA_ENABLE_OPENMP})
   if(NOT (${ENABLE_OPENMP} OR ${ENABLE_CUDA}))
           message(STATUS "ENABLE_OPENMP is set to " ${ENABLE_OPENMP} " ENABLE_CUDA " ${ENABLE_CUDA})
	   message(STATUS "The TPLs have been built without enabling one of the OPENMP or CUDA feature required for USE_RAJA option: please rebuild the RAJA library ?")
   endif()
   
   #if(DEFINED CAMP_DIR)
   #    find_package(camp REQUIRED
   #	       PATHS ${CAMP_DIR}
   #                NO_DEFAULT_PATH)
   #else()
   #   # CAMP is required when using RAJA. Consider using find_package when RAJA_DIR not defined, that is GUIX_INSTALLED_TPLS=ON 
   #   message(STATUS "ProxyApp requires CAMP: CAMP_DIR not set.") 
   #   find_package(camp REQUIRED)
   #endif()

   if(DEFINED RAJA_DIR)
       find_package(RAJA REQUIRED
                    PATHS ${RAJA_DIR}
                    NO_DEFAULT_PATH)
       find_package(camp REQUIRED
	       PATHS ${CAMP_DIR}
	       	     ${camp_DIR}
                    NO_DEFAULT_PATH)
   else()
       message(STATUS "ProxyApp requires RAJA: RAJA_DIR not set.") 
       find_package(camp REQUIRED)
       find_package(RAJA REQUIRED)
   endif()
   get_target_property(CAMP_INCLUDE_DIRS camp INTERFACE_INCLUDE_DIRECTORIES)
   set_target_properties(camp PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CAMP_INCLUDE_DIRS}")
   message("-- camp_DIR: ${camp_DIR} or ${CAMP_DIR}")
   message("-- CAMP_INCLUDE_DIRS: ${CAMP_INCLUDE_DIRS}")
   set (dependencyList ${dependencyList} camp )
   
   get_target_property(RAJA_INCLUDE_DIRS RAJA INTERFACE_INCLUDE_DIRECTORIES)
   set_target_properties(RAJA PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${RAJA_INCLUDE_DIRS}")
   message(STATUS "USE RAJA in RAJA_DIR = ${RAJA_DIR}")      
   message("-- RAJA_INCLUDE_DIRS: ${RAJA_INCLUDE_DIRS}")

   if( ENABLE_HIP )
     set(RAJA_ENABLE_HIP ON CACHE BOOL "" FORCE)
     set(RAJA_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
     include_directories("${ROCM_ROOT_DIR}/include")
   endif()
   
   ################################
   # Umpire
   ################################
   # Umpire is not mandatory for RAJA
   #if(FALSE)
   if(ENABLE_UMPIRE OR ENABLE_CHAI)
	   if(DEFINED UMPIRE_DIR)
		   find_package(umpire REQUIRED
			   PATHS ${UMPIRE_DIR}
			   NO_DEFAULT_PATH)
	   else()
		   message(STATUS "ProxyApp requires Umpire: UMPIRE_DIR not set for finding package.")
		   find_package(umpire REQUIRED)
	   endif()
	   
	   get_target_property(UMPIRE_INCLUDE_DIRS umpire INTERFACE_INCLUDE_DIRECTORIES)
   	   set_target_properties(umpire PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${UMPIRE_INCLUDE_DIRS}")
	   message(STATUS " USE RAJA: RAJA will be compiled with UMPIRE in UMPIRE_DIR = ${UMPIRE_DIR}")
	   message("-- UMPIRE_INCLUDE_DIRS: ${UMPIRE_INCLUDE_DIRS} / umpire_DIR "${umpire_DIR})
	   set (extraComponentsLinkList ${extraComponentsLinkList} umpire )
   endif()
   #endif()
   
   ################################
   # CHAI
   ################################
   # CHAI is not a requisite for RAJA. Umpire is a requisite for CHAI 
   if(ENABLE_CHAI)
	   set(ENABLE_UMPIRE ON CACHE BOOL "" FORCE)
	   if(DEFINED CHAI_DIR)
		   find_package(chai REQUIRED 
			   PATHS ${CHAI_DIR} 
			   NO_DEFAULT_PATH)
	   else()
		   message(STATUS "ProxyApp requires CHAI: CHAI_DIR not set for finding package.")
		   find_package(chai REQUIRED)
	   endif()
	   get_target_property(CHAI_INCLUDE_DIRS chai INTERFACE_INCLUDE_DIRECTORIES)
	   set_target_properties(chai PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CHAI_INCLUDE_DIRS}")
	   message(STATUS " ENABLE_CHAI in CHAI_DIR = ${CHAI_DIR}")
	   message("-- CHAI_INCLUDE_DIRS: ${CHAI_INCLUDE_DIRS} / chai_DIR ${chai_DIR}")
	   set (extraComponentsLinkList ${extraComponentsLinkList} chai )
   endif()
   if(ENABLE_OPENMP)
	   set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
   endif()
   
   ################################
   # CALIPER
   ################################
if(ENABLE_CALIPER)
	    #find_package(adiak REQUIRED)
    if(EXISTS ${CALIPER_DIR})
      find_package(caliper REQUIRED
                 PATHS ${CALIPER_DIR})
    else()
      message(STATUS "CALIPER_DIR not set for finding package")
      find_package(caliper REQUIRED)
    endif()
    get_target_property(CALIPER_INCLUDE_DIRS caliper INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(caliper PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CALIPER_INCLUDE_DIRS}")
    message(STATUS "Using caliper from ${CALIPER_DIR}")
    set (extraComponentsLinkList ${extraComponentsLinkList} caliper )
else()
    message(STATUS "Not using caliper.")
endif()

if(USE_LVARRAY_EXT)
	# Find package
	if(NOT EXISTS ${LVARRAY_DIR})
		message(FATAL_ERROR "LVARRAY_DIR must be set to a defined directory to find the lvarray lib. LVARRAY_DIR is now set to "${LVARRAY_DIR})
	endif()
	find_package(LVARRAY REQUIRED)
	get_target_property(LVARRAY_INCLUDE_DIRS lvarray INTERFACE_INCLUDE_DIRECTORIES)
	set_target_properties(lvarray PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${LVARRAY_INCLUDE_DIRS}")
	message(STATUS "-- Find package for lvarray - LVARRAY_INCLUDE_DIRS " ${LVARRAY_INCLUDE_DIRS})
else()	
   add_subdirectory( LvArray )
   include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LvArray/src)
endif()
   message(STATUS "USE LvArray in LVARRAY_DIR = ${LVARRAY_DIR}")

   if( ENABLE_HIP )
     set (dependencyList ${dependencyList} blt::hip blt::hip_runtime)
   endif()

   set (dependencyList ${dependencyList} lvarray RAJA)
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
else()
  set (USE_LVARRAY OFF CACHE BOOL "" FORCE)
  set (USE_VECTOR OFF CACHE BOOL "" FORCE)
  
  message ("-- USE_KOKKOS flag is set to be true and KOKOS is now enabled")
  set (SOLVER "Kokkos" CACHE PATH "" FORCE)

  # setup KOKOS
  if(DEFINED KOKKOS_DIR)
	  find_package(Kokkos REQUIRED
                    PATHS ${KOKKOS_DIR}
                    NO_DEFAULT_PATH)
  else()
	  find_package(Kokkos REQUIRED)
  endif()
  message("-- Found the KOKKOS library in ${KOKKOS_DIR}")
  
  get_target_property(KOKKOS_INCLUDE_DIRS Kokkos::kokkos INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(KOKKOS_LINK_LIBS Kokkos::kokkos INTERFACE_LINK_LIBRARIES)
  
  message(STATUS "KOKKOS_INCLUDE_DIRS: ${KOKKOS_INCLUDE_DIRS}")     
  message("-- KOKKOS_LINK_LIBS: ${KOKKOS_LINK_LIBS}")

  #get_target_property(KOKKOSCORE_LINK_OPTION Kokkos::kokkoscore INTERFACE_LINK_OPTIONS)
  #message(STATUS "KOKKOSCORE_LINK_OPTION: ${KOKKOSCORE_LINK_OPTION}")

  set (dependencyList ${dependencyList} Kokkos::kokkos)
  if(ENABLE_OPENMP)
	  set (extraComponentsLinkList ${extraComponentsLinkList} openmp )
  endif()
  if( ENABLE_HIP )
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "" FORCE)
    include_directories("${ROCM_ROOT_DIR}/include")
  endif()
  message (STATUS "KOKKOS DEPENDENCY LIST= ${dependencyList}")
  message (STATUS "--KOKKOS enabled: Kokkos_ENABLE_CUDA = ${Kokkos_ENABLE_CUDA} and Kokkos_ENABLE_OPENMP = ${Kokkos_ENABLE_OPENMP}")

endif()
  # Using OMP
  #if(((USE_KOKKOS OR USE_RAJA) AND ENABLE_OPENMP) OR USE_OMP)
if(ENABLE_OPENMP OR USE_OMP)
	message(STATUS "-- USE_OMP: enforcing the flags for OpenMP - Flag = " ${ENABLE_OPENMP})
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" CACHE STRING "" FORCE)
  
  
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OpenMP_C_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}" CACHE STRING "" FORCE)
  
  set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${OpenMP_C_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}" CACHE STRING "" FORCE)
  
  set (CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${OpenMP_C_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${OpenMP_CXX_FLAGS}" CACHE STRING "" FORCE)
  
  set (CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} ${OpenMP_C_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${OpenMP_CXX_FLAGS}" CACHE STRING "" FORCE)
  
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} ${OpenMP_Fortran_FLAGS}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL} ${OpenMP_Fortran_FLAGS}" CACHE STRING "" FORCE)
endif()
message(STATUS "-- The dependencies for building the proxy-app are: " ${dependencyList} " extraComponentsLinkLists = "${extraComponentsLinkList})
if((USE_KOKKOS OR USE_RAJA) AND (ENABLE_OPENMP OR ENABLE_MPI))
	# SettingFlag associated to the config.. used for the name of the executable
	## What programming on the host
	
	if(ENABLE_OPENMP)
		set(flag_hmodel "OMP")
	elseif(ENABLE_MPI)
		set(flag_hmodel "MPI")
		#else()
		#set(flag_hmodel "SERIAL")
	endif()
	
	set(SOLVER "${SOLVER}_${flag_hmodel}" CACHE PATH "" FORCE)
endif()
# The solver's flag
message(STATUS "The flag related to the enabled model on the host is: ${flag_hmodel} - The solver is: " ${SOLVER})

message(STATUS "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")

configure_file( utils/commonConfig.hpp.in
                ${CMAKE_BINARY_DIR}/src/commonConfig.hpp)

add_subdirectory( sem )
add_subdirectory( fd )
if ( addTests)
  add_subdirectory( tests )
endif()
