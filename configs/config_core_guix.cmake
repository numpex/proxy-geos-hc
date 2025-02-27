# This config_core file is included by the TPLs after all the ENABLED options are set. 
# It is also included from the main side in sequential case, after forcing some options 

## C options
#set(CMAKE_C_COMPILER ${CC_ROOT}/gcc${CC_VERSION} CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=${CPU_TUNE_FLAG} -mtune=${CPU_TUNE_FLAG}" CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

## C++ options
#set(CMAKE_CXX_COMPILER ${CXX_ROOT}/g++${CC_VERSION} CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=${CPU_TUNE_FLAG} -mtune=${CPU_TUNE_FLAG}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

## Fortran options
#set(CMAKE_Fortran_COMPILER ${GFORTRAN_ROOT}/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=${CPU_TUNE_FLAG} -mtune=${CPU_TUNE_FLAG}" CACHE STRING "")
#set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")


## OpenMP options
#set(OpenMP_Fortran_FLAGS "-qsmp=omp" CACHE STRING "")
#set(OpenMP_Fortran_LIB_NAMES "" CACHE STRING "")

## MPI options - possible to use the option ENABLE_FIND_MPI with BLT to detect automatically the options and flags
if ( ENABLE_MPI)
     set(ENABLE_MPI ON CACHE BOOL "")
     set(MPI_C_COMPILER         ${MPI_ROOT}/bin/mpicc  CACHE PATH "")
     set(MPI_CXX_COMPILER       ${MPI_ROOT}/bin/mpicxx CACHE PATH "")
     set(MPI_Fortran_COMPILER   ${MPI_ROOT}/bin/mpifort CACHE PATH "")
     set(MPI_FC                 ${MPI_ROOT}/bin/mpifort CACHE PATH "")
     set(MPIEXEC                ${MPI_ROOT}/bin/mpirun  CACHE STRING "")
     set(MPIEXEC_EXECUTABLE     ${MPI_ROOT}/bin/mpirun  CACHE STRING "")
     set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")
endif()

if(ENABLE_CUDA)
  # Cuda options
  set(ENABLE_CUDA_NVTOOLSEXT ON CACHE BOOL "") #Not used
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "")
    #set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "") 
    set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_COMPUTE} CACHE STRING "")
    
    set(CMAKE_CUDA_STANDARD 17 CACHE STRING "") 
    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "")
    ### The inclusion of -std=c++17 is a workaround for a cuda10/gcc8 bug ###
    set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -Xcompiler -std=c++17" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3 -Xcompiler -march=${CPU_TUNE_FLAG} -Xcompiler -mtune=${CPU_TUNE_FLAG}" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
    set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")

    # Uncomment this line to make nvcc output register usage for each kernel.
    # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --resource-usage" CACHE STRING "" FORCE)
endif()

if(ENABLE_HIP)
  set(CMAKE_CXX_COMPILER "${ROCM_ROOT_DIR}/bin/hipcc" CACHE PATH "")
  set(ROCM_PATH ${ROCM_ROOT_DIR} CACHE PATH "")
  set(HIP_ROOT_DIR "${ROCM_ROOT_DIR}/hip" CACHE PATH "")
  set(HIP_PATH "${ROCM_ROOT_DIR}/bin" CACHE PATH "")
  set(GPU_TARGETS "${COMP_ARCH}" CACHE STRING "")
  message(STATUS "======== ROCM_ROOT_DIR=${ROCM_ROOT_DIR}")
endif()

