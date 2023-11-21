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

if (ENABLE_LVARRAY)
  set( dependencyList lvarray RAJA chai )
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

if (ENABLE_KOKKOS)
  set (kokkosList "-L${Kokkos_ROOT}/../../ -lkokkoscontainers -lkokkoscore -lkokkossimd")
  set( dependencyList ${dependencyList} ${kokkosList} )
endif()

if (ENABLE_CUDA)
  if(ENABLE_KOKKOS)
    message( STATUS "kokkosinclude=${Kokkos_ROOT}/../../../include")
    message( STATUS "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
    blt_add_executable( NAME FDTDKOKKOS.exe 
                SOURCES FDTDKokkos.cpp
		HEADERS ../utils/dataType.hpp ../utils/utils.hpp ./utilsFDTD/FDTDutils.hpp
                DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
    target_include_directories(FDTDKOKKOS.exe PUBLIC ${CMAKE_BINARY_DIR}/src)
    target_include_directories(FDTDKOKKOS.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../utils)
  endif()
  if(ENABLE_LVARRAY)
    blt_add_executable( NAME FDTDRAJA.exe 
		        SOURCES FDTDRAJA.cpp
		        HEADERS ../utils/dataType.hpp ../utils/utils.hpp ./utilsFDTD/FDTDutils.hpp
                        DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
    target_include_directories(FDTDRAJA.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)
    target_include_directories(FDTDRAJA.exe PUBLIC ${CMAKE_BINARY_DIR}/src)
    target_include_directories(FDTDRAJA.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../utils)
  endif()
elseif(ENABLE_KOKKOS AND ENABLE_OPENMP)
    message( STATUS "kokkosinclude=${Kokkos_ROOT}/../../../include")
    blt_add_executable( NAME FDTDKOKKOS.exe 
                        SOURCES FDTDKokkos.cpp
		        HEADERS ../utils/dataType.hpp ../utils/utils.hpp ./utilsFDTD/FDTDutils.hpp
                        DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
    target_include_directories(FDTDKOKKOS.exe PUBLIC "${Kokkos_ROOT}/../../../include")
    target_include_directories(FDTDKOKKOS.exe PUBLIC ${CMAKE_BINARY_DIR}/src)
    target_include_directories(FDTDKOKKOS.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../utils)
    target_include_directories(FDTDKOKKOS.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utilsFDTD)
elseif(ENABLE_OPENMP)
     blt_add_executable( NAME FDTDOMP.exe 
                         SOURCES FDTDOMP.cpp
		         HEADERS ../utils/dataType.hpp ../utils/utils.hpp ./utilsFDTD/FDTDutils.hpp
                         DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
     target_include_directories(FDTDOMP.exe PUBLIC ${CMAKE_BINARY_DIR}/src)
     target_include_directories(FDTDOMP.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../utils)
     target_include_directories(FDTDOMP.exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/utilsFDTD)
endif()

