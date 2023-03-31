
set(HEADERS
    utils.hpp
    commonConfig.hpp
    commonMacro.hpp
    )
include_directories(../mesh)

set(SEM_USE_CALIPER FALSE)

if (DEFINED caliper_DIR)
  set(SEM_USE_CALIPER TRUE)
endif()

configure_file(commonConfig.hpp.in commonConfig.hpp)


