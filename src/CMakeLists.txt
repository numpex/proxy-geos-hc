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
find_package(OpenMP REQUIRED)
target_link_libraries(sem.exe
  PRIVATE
  finiteElements
  mesh
  solver
  OpenMP::OpenMP_CXX
  )
install(TARGETS sem.exe
        RUNTIME DESTINATION bin)