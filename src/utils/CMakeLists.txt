
set(HEADERS
    utils.hpp
    commonConfig.hpp
    commonMacro.hpp
    )
include_directories(../mesh)


set(SEM_USE_CALIPER FALSE)

if(EXISTS $ENV{CALIPER_DIR}/caliper-config.cmake)
  set(caliper_DIR "$ENV{CALIPER_DIR}" CACHE PATH "" FORCE)
  message("-- Found the caliper library in ${caliper_DIR}")
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
  message("-- Caliper_DIR is not the correct path and caliper is not enabled")
endif()


set(SEM_USE_RAJA FALSE)

if(EXISTS $ENV{RAJA_DIR}/raja-config.cmake)
  set(RAJA_DIR "$ENV{RAJA_DIR}" CACHE PATH "" FORCE)
  message("-- Found the raja library in ${RAJA_DIR}")
  set(SEM_USE_RAJA TRUE)
  find_package(RAJA)
  target_link_libraries(sem.exe
    PRIVATE 
    RAJA
    )
elseif(NOT DEFINED RAJA_DIR ) 
  message("-- RAJA_DIR is not defined and RAJA is not enabled")
else ()
  message("-- RAJA_DIR is not the correct path and RAJA is not enabled")
endif()

configure_file(commonConfig.hpp.in commonConfig.hpp)


