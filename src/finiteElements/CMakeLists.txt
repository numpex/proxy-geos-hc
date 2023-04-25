set(SRCS
	QkGL.cpp
    )

set(HEADERS
	QkGL.hpp
    )

include_directories(../utils)
include_directories(${PROJECT_BINARY_DIR}/utils)
include_directories(${caliper_INCLUDE_DIR})

add_library(finiteElements ${SRCS} ${HEADERS})

