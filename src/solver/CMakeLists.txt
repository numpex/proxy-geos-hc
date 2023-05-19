message("The current solver is running with ${SOLVER}")

set(SRCS
    solver_${SOLVER}.cpp
    )

set(HEADERS
    solver.hpp
    )

include_directories(../mesh)
include_directories(../finiteElements)

add_library(solver ${SRCS} ${HEADERS})
