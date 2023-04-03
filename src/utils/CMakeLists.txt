
set(HEADERS
    utils.hpp
    commonConfig.hpp
    commonMacro.hpp
    )
include_directories(../mesh)

set(SEM_USE_CALIPER FALSE)

if(EXISTS ${CALIPER_DIR}/share/cmake/caliper/caliper-config.cmake)
  message("-- Found the caliper library in ${CALIPER_DIR}")
  set(caliper_DIR "$ENV{CALIPER_DIR}/share/cmake/caliper" CACHE PATH "" FORCE)
  set(SEM_USE_CALIPER TRUE)
  find_package(caliper)
  target_link_libraries(sem.exe
    PRIVATE
    caliper
    )
  target_include_directories(sem.exe 
    PRIVATE 
    ${caliper_INCLUDE_DIR}
    )
elseif(NOT DEFINED CALIPER_DIR ) 
  message("-- Caliper_DIR is not defined and caliper is not enabled")
else ()
  message("-- Caliper_DIR is not the current path and caliper is not enabled")
endif()



set(SEM_USE_RAJA FALSE)

if(EXISTS ${RAJA_DIR}/lib/cmake/raja/raja-config.cmake)
  message("-- Found the raja library in ${RAJA_DIR}")
  set(RAJA_DIR "$ENV{RAJA_DIR}" CACHE PATH "" FORCE)
  set(SEM_USE_RAJA TRUE)
  find_package(RAJA
               NO_DEFAULT_PATH
               PATHS ${RAJA_DIR}/lib/cmake/raja)
  target_link_libraries(sem.exe
    PUBLIC
    RAJA
    )
elseif(NOT DEFINED RAJA_DIR ) 
  message("-- RAJA_DIR is not defined and RAJA is not enabled")
else ()
  message("-- RAJA_DIR is not the current path and RAJA is not enabled")
endif()

configure_file(commonConfig.hpp.in commonConfig.hpp)


