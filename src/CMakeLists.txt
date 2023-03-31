cmake_minimum_required(VERSION 3.22)
project(sem)

add_executable(sem.exe main/main.cpp)
add_subdirectory(finiteElements)
add_subdirectory(mesh)
add_subdirectory(solver)
add_subdirectory(utils)

include_directories(finiteElements)
include_directories(mesh)
include_directories(solver)
include_directories(utils)
include_directories(${PROJECT_BINARY_DIR}/utils)

find_package(OpenMP REQUIRED)

if (DEFINED caliper_DIR)
  message("-- Found the caliper library in ${caliper_DIR}")
  find_package(caliper)
  target_link_libraries(sem.exe
    PRIVATE
    caliper
    )
  target_include_directories(sem.exe PRIVATE ${caliper_INCLUDE_DIR})
else ()
  message("-- Caliper is not enabled")
endif()

target_link_libraries(sem.exe
  PRIVATE
  finiteElements
  mesh
  solver
  OpenMP::OpenMP_CXX
  )

install(TARGETS sem.exe
        RUNTIME DESTINATION bin)
