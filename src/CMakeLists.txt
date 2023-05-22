add_executable(sem.exe main/main.cpp)

#-----------------------------------------------------
# Programming Models: OMP and RAJA
#-----------------------------------------------------

# the solver is using sequential mode by default 
set (SOLVER "sequentialVector" CACHE PATH "" FORCE)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "-O3 -std=c++14 -march=native")
set( CMAKE_VERBOSE_MAKEFILE on )

# the solver is using OMP multi-threading if SEM_USE_OMP is ON
if (NOT SEM_USE_OMP)
  message ("-- SEM_USE_OMP flag is set to be false and OMP is not enabled")
else()
  message ("-- SEM_USE_OMP flag is set to be true and OMP is now enabled")
  set (SOLVER "ompVector" CACHE PATH "" FORCE)
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# the solver is using RAJA if SEM_USE_RAJA is ON
if (NOT SEM_USE_RAJA)
  message ("-- SEM_USE_RAJA flag is set to be false and RAJA is not enabled")
elseif (EXISTS $ENV{RAJA_DIR}/raja-config.cmake)
  message ("-- SEM_USE_RAJA flag is set to be true and RAJA is now enabled")
  message("-- Found the RAJA library in $ENV{RAJA_DIR}")
  set (SOLVER "rajaVector" CACHE PATH "" FORCE)
  # RAJA uses OMP
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  # setup RAJA
  find_package(RAJA REQUIRED)
  set (raja_INCLUDE_DIR "$ENV{RAJA_INCLUDE_DIR}")
  include_directories(${raja_INCLUDE_DIR})
  target_include_directories(sem.exe 
    PUBLIC 
    ${raja_INCLUDE_DIR}
  )
 
elseif (NOT DEFINED $ENV{RAJA_DIR}) 
  set (SEM_USE_RAJA OFF CACHE BOOL "")
  message("-- RAJA_DIR is not defined and caliper is not enabled")
else ()
  set (SEM_USE_RAJA OFF CACHE BOOL "")
  message("-- RAJA_DIR is not the correct path and caliper is not enabled")
endif()


#-----------------------------------------------------
# Profiler: CALIPER
#-----------------------------------------------------

# for profiling using caliper
if (NOT SEM_USE_CALIPER)
  message ("-- SEM_USE_CALIPER flag is set to be false and CALIPER is not enabled")
elseif (EXISTS $ENV{CALIPER_DIR}/caliper-config.cmake)
  message ("-- SEM_USE_CALIPER flag is set to be true and CALIPER is now enabled")
  set (caliper_DIR "$ENV{CALIPER_DIR}" CACHE PATH "" FORCE)
  message ("-- Found the caliper library in ${caliper_DIR}")
  find_package(caliper REQUIRED)
  target_link_libraries(sem.exe
    PUBLIC
    caliper
    )
  target_include_directories(sem.exe 
    PUBLIC 
    ${caliper_INCLUDE_DIR}
    )
  include_directories(${caliper_INCLUDE_DIR})
elseif (NOT DEFINED ENV{CALIPER_DIR}) 
  set (SEM_USE_CALIPER FALSE CACHE BOOL "")
  message("-- Caliper_DIR is not defined and caliper is not enabled")
else ()
  set (SEM_USE_CALIPER FALSE CACHE BOOL "")
  message("-- Caliper_DIR is not the correct path and caliper is not enabled")
endif()

include_directories(utils)
include_directories(finiteElements)
include_directories(mesh)
include_directories(solver/${SOLVER})
include_directories(${PROJECT_BINARY_DIR}/src/utils)

add_subdirectory(utils)
add_subdirectory(finiteElements)
add_subdirectory(mesh)
add_subdirectory(solver)

target_link_libraries(sem.exe
  PUBLIC
  finiteElements
  mesh
  solver
  )
if( SEM_USE_RAJA)
  target_link_libraries(sem.exe
  PUBLIC 
  RAJA)
endif()

install(TARGETS sem.exe
        RUNTIME DESTINATION bin)


