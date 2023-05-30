
set(HEADERS
    utils.hpp
    dataType.hpp
    commonConfig.hpp
    commonMacro.hpp
    )
include_directories(../mesh)
configure_file(commonConfig.hpp.in commonConfig.hpp)
