message("-- The Current Solver is Using: ${SOLVER}")

set(SRCS
    solverBase.cpp
    solver${SOLVER}.cpp
    )

set(HEADERS
    solverBase.hpp
    solver${SOLVER}.hpp
    )

include_directories(../mesh)
include_directories(../finiteElements)

add_library(solver ${SRCS} ${HEADERS})
