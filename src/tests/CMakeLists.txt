################################
# Add all the tests
################################

# the followings need "-DUSE_KOKKOS=ON -DENABLE_HIP=ON"
#set (testcase  kokkos_hip)

# the followings need "-DUSE_KOKKOS=ON -DENABLE_CUDA=ON"
#set (testcase  Kokkos)
#set (testcase  Kokkos0)
set (testcase  Kokkos1)
#set (testcase  FDTDKokkos)

# the followings need "-DUSE_RAJA=ON -DENABLE_HIP=ON"
#set (testcase  raja_hip_daxpy)

# the followings need "-DUSE_RAJA=ON -DENABLE_CUDA=ON"
#set (testcase  LvArray)
#set (testcase  LvArrayRAJACudaStencil)
#set (testcase  LvArrayRAJACuda)
#set (testcase  LvArray_cuda)
#set (testcase  RAJAInline)
#set (testcase  Raja_cuda)
#set (testcase  FDTDRAJA)
#set (testcase  raja_cuda_daxpy)

# the followings need "-DENABLE_CUDA=ON"
#set (testcase  cuda)
#set (testcase  cuda_runtime)

# the followings need "-DOMP=ON"
if (USE_OMP) 
  
  LIST(APPEND TESTCASES "omp")
#  set (testcase  omp)
#  set (testcase  FDTDOMP)
  message(STATUS "Building tests with OpenMP: ${TESTCASES}")
endif()

FOREACH (testcase ${TESTCASES})
  set (TEST_EXECUTABLE "test_${testcase}.exe" CACHE PATH "" FORCE)
  
  message ( "\n---------------------------------------------------------------------------------" )
  message ( "-- The Current executable to be built is: ${TEST_EXECUTABLE}")
  message ( "-- CMAKE_SOURCE_DIR : ${CMAKE_SOURCE_DIR}")
  message ( "-- PROJECT_BINARY_DIR :  ${PROJECT_BINARY_DIR}")
  message ( "-- extraComponentsLinkList:  ${extraComponentsLinkList}")
  message ( "-- dependencyList:  ${dependencyList}")
  message ( "---------------------------------------------------------------------------------\n" )
  
  blt_add_executable( NAME ${TEST_EXECUTABLE}
                      SOURCES test_${testcase}.cpp
                      DEPENDS_ON ${extraComponentsLinkList} ${dependencyList})
  
  install(TARGETS ${TEST_EXECUTABLE}
          RUNTIME DESTINATION bin)
  add_test(NAME ${TEST_EXECUTABLE}
           COMMAND ${TEST_EXECUTABLE})
ENDFOREACH()
