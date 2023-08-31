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

blt_add_executable( NAME test_LvArray.exe 
                    SOURCES test_LvArray.cpp
                    DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})

target_include_directories(test_LvArray.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)

if (ENABLE_CUDA)
  if(ENABLE_KOKKOS)
	  message( STATUS "kokkosinclude=${Kokkos_ROOT}/../../../include")
    blt_add_executable( NAME test_KOKKOS.exe 
                        SOURCES test_Kokkos.cpp
                        DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
    blt_add_executable( NAME test1_KOKKOS.exe 
                        SOURCES test1_Kokkos.cpp
                        DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
    target_include_directories(test_KOKKOS.exe PUBLIC "${Kokkos_ROOT}/../../../include")
  endif()
  if(ENABLE_LVARRAY)
    blt_add_executable( NAME test_LvArrayRAJACuda.exe 
                        SOURCES test_LvArrayRAJACuda.cpp
                        DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
		target_include_directories(test_LvArrayRAJACuda.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)
    blt_add_executable( NAME    test_cuda.exe
                        SOURCES test_cuda.cpp
                        DEPENDS_ON cuda)
    target_include_directories(test_cuda.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)

    blt_add_executable( NAME    test_cuda_runtime.exe
                        SOURCES test_cuda_runtime.cpp
                        DEPENDS_ON cuda_runtime)
    target_include_directories(test_cuda_runtime.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)

    blt_add_executable( NAME    test_Raja_cuda.exe
	                      SOURCES test_Raja_cuda.cpp
                        DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
    target_include_directories(test_Raja_cuda.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)

    blt_add_executable( NAME    test_LvArray_cuda.exe
	                      SOURCES test_LvArray_cuda.cpp
                        DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
    target_include_directories(test_LvArray_cuda.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)
  endif()
else()
  if(ENABLE_LVARRAY)
    blt_add_executable( NAME test_LVARRAY.exe 
                        SOURCES test_LvArray.cpp
                        DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
    target_include_directories(test_LVARRAY.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)
  endif()
  if(ENABLE_KOKKOS)
    message( STATUS "kokkosinclude=${Kokkos_ROOT}/../../../include")
    blt_add_executable( NAME test_KOKKOS.exe 
                        SOURCES test_Kokkos.cpp
                        DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
    target_include_directories(test_KOKKOS.exe PUBLIC "${Kokkos_ROOT}/../../../include")
  endif()
endif()

