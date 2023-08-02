################################
# Add all the tests
################################
set( extraComponentsLinkList "" )

if( ENABLE_OPENMP )
  set( extraComponentsLinkList ${extraComponentsLinkList} openmp )
endif()

if( ENABLE_CUDA )
  set( extraComponentsLinkList ${extraComponentsLinkList} cuda )
endif()

set( dependencyList lvarray RAJA chai )

if( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper )

  if (TARGET adiak::adiak)
    set( dependencyList ${dependencyList} adiak::adiak )
  endif()

  if (TARGET adiak)
    set( dependencyList ${dependencyList} adiak )
  endif()
endif()

#-------------------------------------------------------
# SEM proxy application setup
#-------------------------------------------------------

set (sem_headers
     utils/commonMacro.hpp
     utils/dataType.hpp
     utils/utils.hpp
     finiteElements/QkGL.hpp
     mesh/simpleMesh.hpp
     solver/solverBase.hpp
     solver/solverRaja.hpp
     main/semproxy.hpp)

set (sem_sources
     finiteElements/QkGL.cpp
     mesh/simpleMesh.cpp
     solver/solverBase.cpp
     solver/solverRaja.cpp
     main/semproxy.cpp
     main/main.cpp)



configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/utils/commonConfig.hpp.in
                ${CMAKE_BINARY_DIR}/include/commonConfig.hpp)

message("--------------------------------------------------CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
message("--------------------------------------------------CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")

blt_add_executable( NAME sem.exe 
                    SOURCES ${sem_sources}
		    HEADERS ${sem_headers}
		    DEPENDS_ON ${extraComponentsLinkList} ${dependencyList} )
target_include_directories(sem.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../LvArray/src)
target_include_directories(sem.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utils)
target_include_directories(sem.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/finiteElements)
target_include_directories(sem.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mesh)
target_include_directories(sem.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/solver)
target_include_directories(sem.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/main)

