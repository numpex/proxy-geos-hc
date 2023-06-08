 #include "../src/mesh/simpleMesh.hpp"
 #include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(pyMesh, m) 
{
    // optional module docstring
    m.doc() = "pybind11 example plugin";

    // bindings to simpleMesh class
    py::class_<simpleMesh>(m, "simpleMesh")
        .def(py::init<int , int , float , float  , int  >())
        .def("getNx", &simpleMesh::getNx);
}