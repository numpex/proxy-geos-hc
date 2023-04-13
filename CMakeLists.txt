cmake_minimum_required(VERSION 3.8)

project(sem)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3 -std=c++14 " CACHE STRING "")

set(SEM_USE_CALIPER TRUE CACHE BOOL "")
set(SEM_USE_RAJA TRUE CACHE BOOL "")
set(SEM_USE_OMP TRUE CACHE BOOL "")

set (SOLVER "sequentialVector" CACHE PATH "" FORCE)
if (SEM_USE_OMP)
  set (SOLVER "ompVector" CACHE PATH "" FORCE)
endif()

add_subdirectory(src)


