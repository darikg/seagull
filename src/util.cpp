#include "util.hpp"

std::vector<Point3> array_to_points_3(const py::array_t<double> &verts) {
    auto v = verts.unchecked<2>();
    if (v.shape(1) != 3) {
        throw std::runtime_error("vertices need to be 3 dimensional");
    }
    const ssize_t nv = v.shape(0);
    std::vector<Point3> points;
    points.reserve(nv);
    for (ssize_t i = 0; i < nv; i++) {
        points.emplace_back(Point3(v(i, 0), v(i, 1), v(i, 2)));
    }
    return points;
}

std::vector<Point2> array_to_points_2(const py::array_t<double> &verts) {
    auto v = verts.unchecked<2>();
    if (v.shape(1) != 2) {
        throw std::runtime_error("vertices need to be 2 dimensional");
    }
    const ssize_t nv = v.shape(0);
    std::vector<Point2> points;
    points.reserve(nv);
    for (ssize_t i = 0; i < nv; i++) {
        points.emplace_back(Point2(v(i, 0), v(i, 1)));
    }
    return points;
}
