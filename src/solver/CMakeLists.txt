cmake_minimum_required(VERSION 3.22)
project(sem)

add_subdirectory(${SOLVER})
include_directories(${SOLVER})
