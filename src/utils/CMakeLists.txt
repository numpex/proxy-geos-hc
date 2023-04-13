
set(HEADERS
    utils.hpp
    commonConfig.hpp
    commonMacro.hpp
    )
include_directories(../mesh)
message("SEM_USE_CALIPER is ${SEM_USE_CALIPER}")
configure_file(commonConfig.hpp.in commonConfig.hpp)
