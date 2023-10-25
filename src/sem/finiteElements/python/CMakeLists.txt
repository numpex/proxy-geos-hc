set (SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/finiteElements/QkGL.cpp")
set (HEADER_FILES "${CMAKE_SOURCE_DIR}/src/finiteElements/QkGL.hpp")
set (PYTHON_FILES "${CMAKE_SOURCE_DIR}/src/finiteElements/python/QkGL.cpp" )

pybind11_add_module(pyQkGL
	                  ${SOURCE_FILES}
	                  ${HEADER_FILES}
	                  ${PYTHON_FILES}
                   )
target_link_libraries(pyQkGL PUBLIC)

install(TARGETS pyQkGL
        LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
       )
