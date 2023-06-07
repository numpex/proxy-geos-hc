 #include "../src/mesh/simpleMesh.hpp"
 #include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(pyMesh, m) 
{
    // optional module docstring
    m.doc() = "pybind11 example plugin";

    // bindings to simpleMesh class
    py::class_<simpleMesh>(m, "simpleMesh")
        .def(py::init<const int & , const int &, const float & , const float & , const int &  >())
        .def("getNumberOfNodes",&simpleMesh::getNumberOfNodes)
        .def("getNx", &simpleMesh::getNx, py::return_value_policy::copy);
}