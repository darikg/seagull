#pragma once
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>


namespace py = pybind11;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::Point_2                 Point2;
typedef Kernel::Point_3                 Point3;
typedef Kernel::Vector_2                Vector2;
typedef Kernel::Vector_3                Vector3;

typedef CGAL::Surface_mesh<Point3>      Mesh3;
typedef Mesh3::Vertex_index             V3;
typedef Mesh3::Face_index               F3;
typedef Mesh3::Halfedge_index           H3;
typedef Mesh3::Edge_index               E3;

typedef CGAL::Surface_mesh<Point2>      Mesh2;
typedef Mesh2::Vertex_index             V2;
typedef Mesh2::Face_index               F2;
typedef Mesh2::Halfedge_index           H2;
typedef Mesh2::Edge_index               E2;
