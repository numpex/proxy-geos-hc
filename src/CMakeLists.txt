#---------------------------------------------------------------
# The Solver Implementation is with Programming Models:
#
#  "OMP" and "RAJA"
#
#  the current implementations include:
#
#  (0) the default is using sequential mode
#  (1) if SEM_USE_OMP is ON it is using OMP to parallelise loops
#  (2) if SEM_USE_RAJA is ON it is using RAJA with vectors
#  (3) if SEM_USE_LVARRAY is ON it is using RAJA with LVARRAY
#
#---------------------------------------------------------------

add_executable(sem.exe main/main.cpp)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "-O3 -std=c++17  " CACHE STRING "" FORCE)

#---------------------------------------------------------------
# (0) the solver is using sequential mode by default
#---------------------------------------------------------------

# the solver is using sequential mode by default 
set (SOLVER "Sequential" CACHE PATH "" FORCE)


#---------------------------------------------------------------
# (1) the solver is using OMP multi-threading if SEM_USE_OMP is ON
#---------------------------------------------------------------

if (NOT SEM_USE_OMP)
  message ("-- SEM_USE_OMP flag is set to be false and OMP is not enabled")
else()
  message ("-- SEM_USE_OMP flag is set to be true and OMP is now enabled")
  set (SOLVER "OMP" CACHE PATH "" FORCE)
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



#-------------------------------------------------------------------------------------
# (2) the solver is using RAJA with vectors if SEM_USE_RAJA and SEM_USE_VECTOR are ON
#-------------------------------------------------------------------------------------

# the solver is using RAJA if SEM_USE_RAJA is ON
if (NOT SEM_USE_RAJA)
  message ("-- SEM_USE_RAJA flag is set to be false and RAJA is not enabled")
elseif (EXISTS $ENV{RAJA_DIR}/raja-config.cmake)
  message ("-- SEM_USE_RAJA flag is set to be true and RAJA is now enabled")
  message("-- Found the RAJA library in $ENV{RAJA_DIR}")
  set (SOLVER "Raja" CACHE PATH "" FORCE)

  # RAJA uses OMP
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

  # setup RAJA
  find_package(RAJA REQUIRED)
  set (raja_INCLUDE_DIR "${RAJA_DIR}/../../../include")

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

#-------------------------------------------------------------------------------------
# (3) the solver is using KOKOS with vectors if SEM_USE_KOKKOS and SEM_USE_VECTOR are ON
#-------------------------------------------------------------------------------------

# the solver is using KOKKOS if SEM_USE_KOKKOS is ON
if (NOT SEM_USE_KOKKOS)
  message ("-- SEM_USE_KOKKOS flag is set to be false and KOKKOS is not enabled")
elseif (EXISTS $ENV{KOKKOS_DIR}/KokkosConfig.cmake)
  message ("-- SEM_USE_KOKKOS flag is set to be true and KOKOS is now enabled")
  message("-- Found the KOKKOS library in $ENV{KOKKOS_DIR}")
  set (SOLVER "Kokkos" CACHE PATH "" FORCE)

  # KOKOS uses OMP
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

  # setup KOKOS
  set(Kokkos_ROOT "$ENV{KOKKOS_DIR}")
  find_package(Kokkos REQUIRED)
  set (Kokkos_INCLUDE_DIRS "${Kokkos_ROOT}/../../../include")
  include_directories(${Kokkos_INCLUDE_DIRS})
  target_include_directories(sem.exe 
    PUBLIC 
    ${Kokkos_INCLUDE_DIRS}
  )

elseif (NOT DEFINED $ENV{KOKKOS_DIR}) 
  set (SEM_USE_KOKKOS OFF CACHE BOOL "")
  message("-- KOKKOS_DIR is not defined and caliper is not enabled")
else ()
  set (SEM_USE_KOKKOS OFF CACHE BOOL "")
  message("-- KOKKOS_DIR is not the correct path and caliper is not enabled")
endif()


#-------------------------------------------------------------------------------------
# (3) the solver is using RAJA with LVARRAY if SEM_USE_RAJA and SEM_USE_LVARRAY are ON
#-------------------------------------------------------------------------------------

if (NOT SEM_USE_LVARRAY)
  message ("-- SEM_USE_LVARRAY flag is set to be false and RAJA is not enabled")
elseif (EXISTS $ENV{LVARRAY_DIR}/lvarray-config.cmake)
  message ("-- SEM_USE_LVARRAY flag is set to be true and LVARRAY is now enabled")
  message("-- Found the LVARRAY library in $ENV{LVARRAY_DIR}")

  # LvArray uses OMP
  find_package(OpenMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

  # LvArray uses  RAJA
  find_package(RAJA REQUIRED)
  set (raja_INCLUDE_DIR "${RAJA_DIR}/../../../include")
  include_directories(${raja_INCLUDE_DIR})
  target_include_directories(sem.exe  
                            PUBLIC
                            ${raja_INCLUDE_DIR}
                            )

  # LvArray uses umpire
  find_package(UMPIRE REQUIRED)
  set (umpire_INCLUDE_DIR "$ENV{UMPIRE_INCLUDE_DIR}")
  include_directories(${umpire_INCLUDE_DIR})
  target_include_directories(sem.exe 
                            PUBLIC 
                            ${umpire_INCLUDE_DIR}
                            )
  # LvArray uses camp
  find_package(CAMP REQUIRED)
  set (camp_INCLUDE_DIR "$ENV{CAMP_INCLUDE_DIR}")
  include_directories(${camp_INCLUDE_DIR})
  target_include_directories(sem.exe 
                            PUBLIC 
                            ${camp_INCLUDE_DIR}
                            )

  # LvArray uses chai
  #find_package(CHAI REQUIRED)
  #set (chai_INCLUDE_DIR "$ENV{CHAI_INCLUDE_DIR}")
  #include_directories(${chai_INCLUDE_DIR})
  #target_include_directories(sem.exe 
  #                          PUBLIC 
  #                          ${chai_INCLUDE_DIR}
  #                          )


  # LvArray 
  find_package(LVARRAY REQUIRED)
  set (lvarray_INCLUDE_DIR "${LVARRAY_DIR}/../../../include")
  include_directories(${lvarray_INCLUDE_DIR})
  target_include_directories(sem.exe 
                            PUBLIC 
                            ${lvarray_INCLUDE_DIR}
                            )

  target_link_libraries(sem.exe
                        PUBLIC 
                        umpire
                        camp
  #                      chai
                        lvarray)

elseif (NOT DEFINED $ENV{LVARRAY_DIR}) 
  set (SEM_USE_LVARRAY OFF CACHE BOOL "")
  message("-- LVARRAY_DIR is not defined and lvarray is not enabled")
else ()
  set (SEM_USE_LVARRAY OFF CACHE BOOL "")
  message("-- LVARRAY_DIR is not the correct path and lvarray is not enabled")
endif()


#-----------------------------------------------------
# CALIPER
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

#------------------------------------------------------
#  pybind11
#------------------------------------------------------
if (SEM_USE_PYTHON)
  set (python_INCLUDE_DIR "/home/j0024947/anaconda3/include/python3.10/")
  find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

  set(pybind11_DIR "$ENV{PYBIND11_DIR}" CACHE PATH "" FORCE)
  message ("-- Found the pybind11 library in ${pybind11_DIR}")
  find_package(pybind11 REQUIRED)
  set (pybind11_INCLUDE_DIR "$ENV{PYBIND11_INCLUDE_DIR}")

  include_directories(${pybind11_INCLUDE_DIR}
                      ${python_INCLUDE_DIR})
else()
  message ("-- SEM_USE_PYTHON flag is set to be false and PYTHON is not enabled")
endif()


message ( "###########  CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}")
message ( "###########  PROJECT_BINARY_DIR ${PROJECT_BINARY_DIR}")
#-------------------------------------------------------
# SEM
#-------------------------------------------------------
include_directories(utils)
include_directories(finiteElements)
include_directories(mesh)
include_directories(solver)
include_directories(main)
include_directories(${PROJECT_BINARY_DIR}/src/utils)

add_subdirectory(utils)
add_subdirectory(finiteElements)
add_subdirectory(mesh)
add_subdirectory(solver)
add_subdirectory(main)

target_link_libraries(sem.exe
  PUBLIC
  semproxy
  finiteElements
  mesh
  solver
  )
  
if( SEM_USE_RAJA)
  target_link_libraries(sem.exe
  PUBLIC 
  RAJA)
endif()

if( SEM_USE_KOKKOS)
  target_link_libraries(sem.exe
  PUBLIC 
  Kokkos::kokkos)
endif()

if( SEM_USE_LVARRAY)
  target_link_libraries(sem.exe
                        PUBLIC 
                        umpire
                        camp
#                       chai
                        lvarray)
endif()

if (SEM_USE_PYTHON)
  set(PYTHON_LIBRARY_DIR "/home/j0024947/anaconda3/lib/python3.10/site-packages")
  set(PYTHON_EXECUTABLE "/home/j0024947/anaconda3/bin/python3.10")

  add_subdirectory(mesh/python)
  add_subdirectory (finiteElements/python)

endif()

install(TARGETS sem.exe
        RUNTIME DESTINATION bin)

