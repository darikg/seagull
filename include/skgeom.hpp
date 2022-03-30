#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <iostream>

namespace py = pybind11;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT                                          SKGEOM_FT;
typedef Kernel::RT                                          SKGEOM_RT;

typedef Kernel::Point_2                                     Point_2;
typedef Kernel::Point_3                                     Point_3;
