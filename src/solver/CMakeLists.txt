message("-- The Current Solver is Using: ${SOLVER}")

set(SRCS
    solverBase.cpp
    solver${SOLVER}.cpp
    )

set(HEADERS
    solverBase.hpp
    solver${SOLVER}.hpp
    )

add_library(solver ${SRCS} ${HEADERS})
