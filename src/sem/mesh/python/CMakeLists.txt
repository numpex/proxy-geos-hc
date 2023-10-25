set (SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/mesh/simpleMesh.cpp"
                  "${CMAKE_SOURCE_DIR}/src/finiteElements/QkGL.cpp"
                  )
set (HEADER_FILES "${CMAKE_SOURCE_DIR}/src/mesh/simpleMesh.hpp"
                  "${CMAKE_SOURCE_DIR}/src/finiteElements/QkGL.hpp")

set (PYTHON_FILES "${CMAKE_SOURCE_DIR}/src/mesh/python/simpleMesh.cpp" )
pybind11_add_module(pyMesh
	                  ${SOURCE_FILES}
	                  ${HEADER_FILES}
	                  ${PYTHON_FILES}
                   )
target_link_libraries(pyMesh PUBLIC)

install(TARGETS pyMesh
        LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
       )
