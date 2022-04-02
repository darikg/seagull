#include "seagullmesh.hpp"

#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Surface_mesh_shortest_path.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

typedef std::vector<F3>         Faces3;
typedef std::vector<V3>         Verts3;
typedef std::vector<Point3>     Points3;
typedef std::vector<Point2>     Points2;

namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<Kernel::FT>    Barycentric_coordinates;
typedef PMP::Face_location<Mesh3, Kernel::FT>       FaceLoc3;

typedef typename CGAL::AABB_face_graph_triangle_primitive<Mesh3>    AABB_primitive3;
typedef typename CGAL::AABB_traits<Kernel, AABB_primitive3>         AABB_traits3;
typedef typename CGAL::AABB_tree<AABB_traits3>                      AABB_Tree3;


template<typename Point>
py::array_t<double, py::array::c_style> points_d_to_array(const std::vector<Point>& points, const size_t ndims) {
    const size_t np = points.size();
    py::array_t<double, py::array::c_style> points_out({np, ndims});
    auto r = points_out.mutable_unchecked<2>();
    for (auto i = 0; i < np; i++) {
        for (auto j = 0; j < ndims; j++) {
            r(i, j) = CGAL::to_double(points[i][j]);
        }
    }
    return points_out;
}

py::array_t<double, py::array::c_style> points_to_array(const Points3& points) {
    return points_d_to_array<Point3>(points, 3);
}

py::array_t<double, py::array::c_style> points_to_array(const Points2& points) {
    return points_d_to_array<Point2>(points, 2);
}

auto construct_points(const Mesh3& mesh, const std::vector<F3>& faces, const py::array_t<double>& bary_coords) {
    size_t nf = faces.size();
    auto rbc = bary_coords.unchecked<2>();
    size_t nb = rbc.shape(0);
    if (nf != nb) {
        throw std::runtime_error("number of faces doesn't match number of points");
    }

    std::vector<Point3> points;
    points.reserve(nf);

    for (auto i = 0; i < nf; i++) {
        Barycentric_coordinates bc = {rbc(i, 0), rbc(i, 1), rbc(i, 2)};
        FaceLoc3 loc = {faces[i], bc};
        auto pt = PMP::construct_point(loc, mesh);
        points.emplace_back(pt);
    }
    return points_to_array(points);
}

auto locate_points(const Mesh3& mesh, const AABB_Tree3& tree, const py::array_t<double>& points) {
    auto rp = points.unchecked<2>();
    size_t np = rp.shape(0);
    std::vector<F3> faces;
    faces.reserve(np);
    py::array_t<double, py::array::c_style> bary_coords({np, size_t(3)});

    auto rbc = bary_coords.mutable_unchecked<2>();

    for (size_t i = 0; i < np; i++) {
        Point3 point = {rp(i, 0), rp(i, 1), rp(i, 2)};
        FaceLoc3 loc = PMP::locate_with_AABB_tree(point, tree, mesh);
        faces.emplace_back(loc.first);

        for (size_t j = 0; j < 3; j++) {
            rbc(i, j) = loc.second[j];
        }
    }

    return std::make_tuple(faces, bary_coords);
}

// struct Point2_to_Point3 {
//     using key_type = V2;
//     using value_type = Point3;
//     using reference = Point3;
//     using category = boost::readable_property_map_tag;

//     const Mesh2& mesh;

//     Point2_to_Point3(const Mesh2 &mesh) : mesh(mesh) {}

//     friend Point_3 get(const Point2_to_Point3 &map, V2 v) {
//         //auto p = map.mesh.point(v);
//         // return {p[0], p[1], 0};
//         return {0, 0, 0};
//     }
// };

void init_locate(py::module &m) {
    py::module sub = m.def_submodule("locate");

    py::class_<AABB_Tree3>(sub, "AABB_Tree3");

    sub.def("construct_points", &construct_points)
        .def("aabb_tree", [](const Mesh3& mesh) {
            AABB_Tree3 tree;
            PMP::build_AABB_tree(mesh, tree);
            return tree;
        })
        .def("locate_points", &locate_points)
        .def("shortest_path", [](
                const Mesh3& mesh,
                const F3 src_face, const std::vector<double>& src_bc,
                const F3 tgt_face, const std::vector<double>& tgt_bc) {

            using SPTraits = CGAL::Surface_mesh_shortest_path_traits<Kernel, Mesh3>;
            using ShortestPath = CGAL::Surface_mesh_shortest_path<SPTraits>;

            Barycentric_coordinates src_bc_ = {src_bc[0], src_bc[1], src_bc[2]};
            Barycentric_coordinates tgt_bc_ = {tgt_bc[0], tgt_bc[1], tgt_bc[2]};

            ShortestPath shortest_path(mesh);
            shortest_path.add_source_point(src_face, src_bc_);
            std::vector<Point3> points;
            shortest_path.shortest_path_points_to_source_points(tgt_face, tgt_bc_, std::back_inserter(points));

            return points_to_array(points);
        })
    ;
}
