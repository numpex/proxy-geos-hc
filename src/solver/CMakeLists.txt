message("The current solver is running with ${SOLVER}")

set(SRCS
    ${SOLVER}_solver.cpp
    )

set(HEADERS
    solver.hpp
    )

include_directories(../mesh)
include_directories(../finiteElements)

add_library(solver ${SRCS} ${HEADERS})
