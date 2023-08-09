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

if(SEM_USE_LVARRAY)
set( dependencyList lvarray  )
endif()

if(SEM_USE_RAJA )
   if (NOT SEM_USE_LVARRAY)
     set( dependencyList RAJA chai )
   endif()
   if(SEM_USE_LVARRAY)
      set( dependencyList lvarray RAJA chai )
   endif()
endif()

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
message("         ")
message( STATUS "-------------------SEM ${SOLVER}_${PRG} Solver settings-----------------------" )
message( STATUS "-----------------------------------------------------------------------" )

set (sem_headers
     utils/commonMacro.hpp
     utils/dataType.hpp
     utils/utils.hpp
     finiteElements/QkGL.hpp
     mesh/simpleMesh.hpp
     solver/solverBase.hpp
     solver/solver${SOLVER}.hpp
     main/semproxy.hpp)

set (sem_sources
     finiteElements/QkGL.cpp
     mesh/simpleMesh.cpp
     solver/solverBase.cpp
     solver/solver${SOLVER}.cpp
     main/semproxy.cpp
     main/main.cpp)



configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/utils/commonConfig.hpp.in
                ${CMAKE_BINARY_DIR}/src/sem/commonConfig.hpp)
# Install the generated header.
install( FILES ${CMAKE_BINARY_DIR}/src/sem/CommonConfig.hpp
         DESTINATION src/sem )
message(STATUS "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")

blt_add_executable( NAME sem${SOLVER}_${PRG}.exe 
                    SOURCES ${sem_sources}
		    HEADERS ${sem_headers}
		    DEPENDS_ON ${extraComponentsLinkList} ${dependencyList} )
if (SEM_USE_LVARRAY)
	target_include_directories(sem${SOLVER}_${PRG}.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../LvArray/src)
	target_include_directories(sem${SOLVER}_${PRG}.exe PUBLIC ${CMAKE_BINARY_DIR}/include)
endif()
target_include_directories(sem${SOLVER}_${PRG}.exe PUBLIC ${CMAKE_BINARY_DIR}/src/sem)
target_include_directories(sem${SOLVER}_${PRG}.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utils)
target_include_directories(sem${SOLVER}_${PRG}.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/finiteElements)
target_include_directories(sem${SOLVER}_${PRG}.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mesh)
target_include_directories(sem${SOLVER}_${PRG}.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/solver)
target_include_directories(sem${SOLVER}_${PRG}.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/main)

message("         ")
message( STATUS "-------------------END SEM ${SOLVER}_${PRG} Solver settings-------------------" )
message( STATUS "-----------------------------------------------------------------------" )
