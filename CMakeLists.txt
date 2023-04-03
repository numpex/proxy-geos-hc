cmake_minimum_required(VERSION 3.8)

project(sem)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3 -std=c++14 " CACHE STRING "")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CALIPER_DIR "$ENV{CALIPER_DIR}" CACHE PATH "" FORCE)

set(RAJA_DIR "$ENV{RAJA_DIR}" CACHE PATH "" FORCE)

add_subdirectory(src)
