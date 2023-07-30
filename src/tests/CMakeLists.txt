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

blt_add_executable( NAME test_LvArray.exe 
                    SOURCES test_LvArray.cpp
                    DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})

target_include_directories(test_LvArray.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)

if (ENABLE_CUDA)
  blt_add_executable( NAME    test_cuda.exe
                      SOURCES test_cuda.cpp
                      DEPENDS_ON cuda)

  blt_add_executable( NAME    test_cuda_runtime.exe
                      SOURCES test_cuda_runtime.cpp
                      DEPENDS_ON cuda_runtime)

  blt_add_executable( NAME    test_Raja_cuda.exe
	              SOURCES test_Raja_cuda.cpp
                      DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})

  blt_add_executable( NAME    test_LvArray_cuda.exe
	              SOURCES test_LvArray_cuda.cpp
                      DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})

  target_include_directories(test_LvArray_cuda.exe PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src)

endif()

