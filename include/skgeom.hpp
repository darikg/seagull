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

typedef Kernel::Point_2                             Point_2;
typedef Kernel::Point_3                             Point_3;

typedef CGAL::Surface_mesh<Point_3>                 Mesh3;
typedef Mesh3::Vertex_index                         V3;
typedef Mesh3::Face_index                           F3;
typedef Mesh3::Halfedge_index                       H3;
typedef Mesh3::Edge_index                           E3;

typedef CGAL::Surface_mesh<Point_2>                 Mesh2;
typedef Mesh2::Vertex_index                         V2;
typedef Mesh2::Face_index                           F2;
typedef Mesh2::Halfedge_index                       H2;
typedef Mesh2::Edge_index                           E2;


py::array_t<double, py::array::c_style> points_to_array(const std::vector<Point_3>& points) {
    // convert points to arrays
    const size_t np = points.size();
    py::array_t<double, py::array::c_style> points_out({np, size_t(3)});
    auto r = points_out.mutable_unchecked<2>();
    for (auto i = 0; i < np; i++) {
        for (auto j = 0; j < 3; j++) {
            r(i, j) = CGAL::to_double(points[i][j]);
        }
    }
    return points_out;
}

std::vector<Point_3> array_to_points_3(const py::array_t<double> &verts) {
    auto v = verts.unchecked<2>();
    if (v.shape(1) != 3) {
        throw std::runtime_error("vertices need to be 3 dimensional");
    }
    const ssize_t nv = v.shape(0);
    std::vector<Point_3> points;
    points.reserve(nv);
    for (ssize_t i = 0; i < nv; i++) {
        points.emplace_back(Point_3(v(i, 0), v(i, 1), v(i, 2)));
    }
    return points;
}

std::vector<Point_2> array_to_points_2(const py::array_t<double> &verts) {
    auto v = verts.unchecked<2>();
    if (v.shape(1) != 2) {
        throw std::runtime_error("vertices need to be 2 dimensional");
    }
    const ssize_t nv = v.shape(0);
    std::vector<Point_2> points;
    points.reserve(nv);
    for (ssize_t i = 0; i < nv; i++) {
        points.emplace_back(Point_2(v(i, 0), v(i, 1)));
    }
    return points;
}

std::vector<Point_3> array_to_points(const Mesh3& mesh, const py::array_t<double> &verts) {
    return array_to_points_3(verts);
}

std::vector<Point_2> array_to_points(const Mesh2& mesh, const py::array_t<double> &verts) {
    return array_to_points_2(verts);
}


constexpr size_t ndims(const Mesh3& mesh) { return 3; }
constexpr size_t ndims(const Mesh2& mesh) { return 2; }