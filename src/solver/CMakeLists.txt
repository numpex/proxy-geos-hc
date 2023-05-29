message("The current solver is running with ${SOLVER}")

set(SRCS
    solver${SOLVER}.cpp
    )

set(HEADERS
    solver.hpp
    )

include_directories(../mesh)
include_directories(../finiteElements)
include_directories(../utils)


add_library(solver ${SRCS} ${HEADERS})
