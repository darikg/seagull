// Python CGAL bindings 
// Author: Wolf Vollprecht <w.vollprecht@gmail.com>

#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_mesh(py::module&);
// void init_properties(py::module&);
// void init_corefine(py::module&);
// void init_meshing(py::module&);
void init_locate(py::module&);

PYBIND11_MODULE(_skgeom, m) {
    m.doc() = "";
    init_mesh(m);
    // init_properties(m);
    // init_corefine(m);
    // init_meshing(m);
    init_locate(m);
}
