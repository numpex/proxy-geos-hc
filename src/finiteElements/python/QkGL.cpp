
 #include "../../finiteElements/QkGL.hpp"
 #include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE( pyQkGL, m )
{
  // optional module docstring
  m.doc() = "QkGL python interface";

  // bindings to simpleMesh class
  py::class_< QkGL >( m, "QkGL" )
    .def( py::init< >());

}
