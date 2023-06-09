 #include "../src/mesh/simpleMesh.hpp"
 #include "../src/finiteElements/QkGL.hpp"
 #include "../src/solver/solver.hpp"
 #include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(pyMesh, m) 
{
    // optional module docstring
    m.doc() = "simpleMesh python interface";

    // bindings to simpleMesh class
    py::class_<simpleMesh>(m, "simpleMesh")
        .def(py::init<int , int , float , float  , int  >())
        .def("getNumberOfNodes",&simpleMesh::getNumberOfNodes)
        .def("getNumberOfElements",&simpleMesh::getNumberOfElements)
        .def("getNx", &simpleMesh::getNx)
        .def("getNy", &simpleMesh::getNy)
        .def("getDx", &simpleMesh::getDx)
        .def("getDy", &simpleMesh::getDy)
        .def("getNumberOfInteriorElements", &simpleMesh::getNumberOfInteriorElements)
        .def("getNumberOfInteriorNodes",&simpleMesh::getNumberOfInteriorNodes)
        .def("getNumberOfBoundaryFaces",&simpleMesh::getNumberOfBoundaryFaces)
        .def("getNumberOfBoundaryNodes",&simpleMesh::getNumberOfBoundaryNodes);
    
    // bindings to simpleMesh class
    py::class_<solver>(m, "solver")
        .def(py::init< >())
        .def("computeFEInit",&solver::computeFEInit);
}