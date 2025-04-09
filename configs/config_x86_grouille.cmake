# SET CONFIG PARAMETERS

get_filename_component(MY_CONFIG_FILE_NAME ${CMAKE_CURRENT_LIST_FILE} NAME)
set(CONFIG_NAME ${MY_CONFIG_FILE_NAME} CACHE PATH "")

site_name(HOST_NAME)

message(STATUS "The CONFIG_NAME is ${CONFIG_NAME} and HOST_NAME = ${HOST_NAME}")
####################################################
# Part 1: Specify below which programming models are enabled
####################################################

## On the Host
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)
## On the device
set(ENABLE_CUDA ON CACHE BOOL "" FORCE)
## Setting ENABLE_OPENMP=ON and ENABLE_CUDA=ON is mandatory for the option USE_RAJA on the proxy-app side


set(ENABLE_MPI OFF CACHE BOOL "" FORCE)
### May need to consider setting ENABLE_FIND_MPI to ON to detect automatically with BLT
set(ENABLE_FIND_MPI OFF CACHE BOOL "" FORCE)

############################################################
# Part 2: Config for the platform and the GCC/G++ or MPI.. compilers
############################################################

## Compute option - GCC/G++ and GFORTRAN compilers
set(CPU_TUNE_FLAG "znver2" CACHE STRING " Flag specifying the mtune and cpu flag for GCC/G++ compiler" FORCE)
set(CC_VERSION "" CACHE STRING "Specify any specfic suffix/version for the GCC/G++ compiler" FORCE)

set(CC_ROOT "/usr/bin" CACHE PATH "Path for GCC compiler" FORCE)
set(CXX_ROOT "/usr/bin" CACHE PATH "Path for G++ compiler" FORCE)
set(GFORTRAN_ROOT "/usr/bin" CACHE PATH "Path for Fortran compiler" FORCE)

## Root for MPI compiler
if(ENABLE_MPI)
  set(MPI_ROOT "/usr/bin" CACHE PATH "Path for MPI compiler" FORCE)
endif()
## Variables for CUDA
if(ENABLE_CUDA)
  set(CUDA_ARCH "sm_80" CACHE STRING "Flag specifying the GPU architecture" FORCE)
  set(CUDA_ARCH_COMPUTE "80" CACHE STRING "Flag specifying the GPU compute capabilities" FORCE)
  set(CUDA_ROOT "/usr/local/cuda" CACHE PATH "The path for the cuda toolkit" FORCE)
  set(CUDA_TOOLKIT_ROOT_DIR ${CUDA_ROOT} CACHE PATH "" FORCE)
  set(CMAKE_CUDA_FLAGS "-lineinfo -g -pg -restrict --expt-relaxed-constexpr --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -std=c++17" CACHE STRING "")
  set(DEVICE "A100" CACHE PATH "" FORCE)
endif()
# Variables for HIP:
if(ENABLE_HIP)
  set(COMP_VER "6.0.0" CACHE STRING "")
  set(COMP_ARCH "gfx908" CACHE STRING "")
  set(ROCM_ROOT_DIR "/opt/rocm-${COMP_VER}" CACHE PATH "")
  set(USE_SEM_INLINE ON CACHE BOOL "" FORCE) #clang doesn't support -fgpu-rdc yet
endif()


###########################################################
# Part 3: Tools for profiling and performance analysis
###########################################################

set(ENABLE_CALIPER OFF CACHE BOOL "")
set(ENABLE_OPENMPTARGET OFF CACHE BOOL "" FORCE)

###########################################################
# Part 4: Include the core config file after all the option are set
###########################################################

# The directory with the config. files
#set(CONFIG_DIR $ENV{proxy_config_root}/configs)

if(EXISTS ${CONFIG_DIR}/config_core.cmake)
        include(${CONFIG_DIR}/config_core.cmake)
else()
        message(FATAL_ERROR "The config_core.cmake is not found in the the config files directory: " ${CONFIG_DIR})
endif()
