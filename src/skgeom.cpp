// Python CGAL bindings 
// Author: Wolf Vollprecht <w.vollprecht@gmail.com>

#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_mesh(py::module&);
void init_properties(py::module&);


PYBIND11_MODULE(_skgeom, m) {
    m.doc() = "";
    init_mesh(m);
    init_properties(m);
}
