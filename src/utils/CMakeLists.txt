
set(HEADERS
    utils.hpp
    commonConfig.hpp
    commonMacro.hpp
    )
include_directories(../mesh)

set(SEM_USE_CALIPER FALSE)

if(EXISTS ${CALIPER_DIR}/caliper-config.cmake)
  message("-- Found the caliper library in ${caliper_DIR}")
  set(caliper_DIR "$ENV{CALIPER_DIR}" CACHE PATH "" FORCE)
  set(SEM_USE_CALIPER TRUE)
  find_package(caliper)
  target_link_libraries(sem.exe
    PRIVATE
    caliper
    )
  target_include_directories(sem.exe PRIVATE ${caliper_INCLUDE_DIR})
else ()
  message("-- Caliper is not enabled")
endif()

configure_file(commonConfig.hpp.in commonConfig.hpp)


