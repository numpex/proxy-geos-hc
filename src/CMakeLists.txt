cmake_minimum_required(VERSION 3.22)
project(sem)
add_executable(sem.exe main/main.cpp)


if (NOT SEM_USE_CALIPER)
  message ("-- SEM_USE_CALIPER flag is set to be false and CALIPER is not enabled")
elseif(EXISTS $ENV{CALIPER_DIR}/caliper-config.cmake)
  set(caliper_DIR "$ENV{CALIPER_DIR}" CACHE PATH "" FORCE)
  message("-- Found the caliper library in ${caliper_DIR}")
  find_package(caliper)
  target_link_libraries(sem.exe
    PUBLIC
    caliper
    )
  target_include_directories(sem.exe 
    PUBLIC 
    ${caliper_INCLUDE_DIR}
    )
elseif(NOT DEFINED CALIPER_DIR ) 
  message("-- Caliper_DIR is not defined and caliper is not enabled")
else ()
  message("-- Caliper_DIR is not the correct path and caliper is not enabled")
endif()


if (NOT SEM_USE_RAJA)
  message ("-- SEM_USE_RAJA flag is set to be false and RAJA is not enabled")
elseif(EXISTS $ENV{RAJA_DIR}/raja-config.cmake)
  set(RAJA_DIR "$ENV{RAJA_DIR}" CACHE PATH "" FORCE)
  message("-- Found the raja library in ${RAJA_DIR}")
  find_package(RAJA)
  target_link_libraries(sem.exe
    PUBLIC 
    RAJA
    )
elseif(NOT DEFINED RAJA_DIR ) 
  message("-- RAJA_DIR is not defined and RAJA is not enabled")
else ()
  message("-- RAJA_DIR is not the correct path and RAJA is not enabled")
endif()


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_subdirectory(utils)
add_subdirectory(finiteElements)
add_subdirectory(mesh)
add_subdirectory(solver)

include_directories(utils)
include_directories(finiteElements)
include_directories(mesh)
include_directories(solver/${SOLVER})
include_directories(${PROJECT_BINARY_DIR}/utils)

target_link_libraries(sem.exe
  PUBLIC
  finiteElements
  mesh
  solver
  )

install(TARGETS sem.exe
        RUNTIME DESTINATION bin)


