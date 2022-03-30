#include "skgeom.hpp"

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Surface_mesh/IO/OFF.h>

namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<Kernel::FT>    Barycentric_coordinates;
typedef PMP::Face_location<Mesh3, Kernel::FT>       Loc3;

typedef typename CGAL::AABB_face_graph_triangle_primitive<Mesh3>    AABB_primitive3;
typedef typename CGAL::AABB_traits<Kernel, AABB_primitive3>         AABB_traits3;
typedef typename CGAL::AABB_tree<AABB_traits3>                      AABB_Tree3;


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

std::vector<Point_3> array_to_points(const Mesh3& mesh, const py::array_t<double> &verts) {
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

std::vector<Point_2> array_to_points(const Mesh2& mesh, const py::array_t<double> &verts) {
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


constexpr size_t ndims(const Mesh3& mesh) { return 3; }
constexpr size_t ndims(const Mesh2& mesh) { return 2; }


template<typename Mesh, typename Point, typename F>
auto construct_points(const Mesh& mesh, const std::vector<F>& faces, const py::array_t<double>& bary_coords) {
    using Loc = typename PMP::Face_location<Mesh, Kernel::FT>;

    size_t nf = faces.size();
    auto rbc = bary_coords.unchecked<2>();
    size_t nb = rbc.shape(0);
    if (nf != nb) {
        throw std::runtime_error("number of faces doesn't match number of points");
    }

    size_t nd = ndims(mesh);
    py::array_t<double, py::array::c_style> points({nf, nd});
    auto rp = points.mutable_unchecked<2>();

    for (auto i = 0; i < nf; i++) {
        Barycentric_coordinates bc = {rbc(i, 0), rbc(i, 1), rbc(i, 2)};
        Loc loc = {faces[i], bc};
        Point pt = PMP::construct_point(loc, mesh);
        for (auto j = 0; j < nd; j ++) {
            rp(i, j) = CGAL::to_double(pt[j]);
        }
    }

    return points;
}


struct Point2_to_Point3 {
    using key_type = V2;
    using value_type = Point_3;
    using reference = Point_3;
    using category = boost::readable_property_map_tag;

    const Mesh2& mesh;

    Point2_to_Point3(const Mesh2 &mesh) : mesh(mesh) {}

    friend Point_3 get(const Point2_to_Point3 &map, V2 v) {
        //auto p = map.mesh.point(v);
        // return {p[0], p[1], 0};
        return {0, 0, 0};
    }
};


struct VertexPointMapWrapper {
    // Used for tracking which verts get moved during remesh, etc
    using key_type = V3;
    using value_type = Point_3;
    using reference = Point_3&;
    using category = boost::read_write_property_map_tag;

    Mesh3::Property_map<V3, Point_3>& points;
    Mesh3::Property_map<V3, bool>& touched;

    VertexPointMapWrapper(Mesh3::Property_map<V3, Point_3>& p, Mesh3::Property_map<V3, bool>& t) : points(p), touched(t) {}

    friend Point_3 get (const VertexPointMapWrapper& map, V3 v) { return map.points[v]; }
    friend void put (const VertexPointMapWrapper& map, V3 v, const Point_3& point) {
        map.points[v] = point;
        map.touched[v] = true;
    }
};


typedef CGAL::Triple<size_t, size_t, size_t> Triangle;

py::array_t<size_t> constrained_contour_pair_mesh(
    const std::vector<Point_3>& p,
    const std::vector<Point_3>& q,
    const std::vector<size_t>& pidx,
    const std::vector<size_t>& qidx,
    const size_t np0,
    const size_t nq0
) {
    const auto np = p.size(), nq = q.size();
    const auto n = pidx.size();
    size_t nf = 0;
    py::array_t<size_t, py::array::c_style> faces({size_t(np + nq), size_t(3)});
    auto r = faces.mutable_unchecked<2>();

    for (auto i = 0; i < (n - 1); i++) {
        const auto p0 = pidx[i], p1 = pidx[i + 1], q0 = qidx[i], q1 = qidx[i + 1];
        // p1 is always > p0 unless p1 is 0
        const auto npi = (p1 != 0) ? p1 - p0 : np - p0;
        const auto nqi = (q1 != 0) ? q1 - q0 : nq - q0;

        // Construct the border of the polygonal face to be triangulated
        // Patch has (npi + 1) P vertices, and (nqi + 1) Q vertices
        // Note that we iterate in reverse order over the q pts
        std::vector<Point_3> polygon;
        polygon.reserve(npi + nqi);
        for (auto j = 0; j <= npi; j++) {
            polygon.emplace_back(p[(p0 + j) % np]);
        }
        for (auto j = nqi + 1; j-- > 0;) {
            polygon.emplace_back(q[(q0 + j) % nq]);
        }

        std::vector<Triangle> patch;
        patch.reserve(npi + nqi - 2);
        PMP::triangulate_hole_polyline(polygon, std::back_inserter(patch));

        // Translate the local patch back into points indices
        for (auto j = 0; j < patch.size(); j++) {
            const auto a = patch[j].first, b = patch[j].second, c = patch[j].third;
            // The q indices are a little hairy because of the reverse ordering.
            // Let v >= (npi + 1) be an index into one of the Q vertices in the patch
            //      q_patch_idx = (v - (npi + 1))  # account for the (npi + 1) P vertices
            // Because of the reverse ordering of the Q points,
            //      q_pts_idx = q0 + (nqi - q_patch_idx)
            //                = q0 + npi + nqi + 1 - v
            r(nf, 0) = (a <= npi) ? ((p0 + a) % np) + np0 : ((q0 + npi + nqi + 1 - a) % nq) + nq0;
            r(nf, 1) = (b <= npi) ? ((p0 + b) % np) + np0 : ((q0 + npi + nqi + 1 - b) % nq) + nq0;
            r(nf, 2) = (c <= npi) ? ((p0 + c) % np) + np0 : ((q0 + npi + nqi + 1 - c) % nq) + nq0;
            nf++;
        }
    }

    return faces;
}


template<typename Mesh, typename Point, typename V, typename F, typename E, typename H>
auto define_mesh(py::module &m, std::string name) {
    return py::class_<Mesh>(m, name.c_str())
        .def(py::init<>())
        .def(py::init([](py::array_t<double> &verts, std::vector<std::vector<size_t>>& faces, const bool orient) {
            Mesh mesh;
            std::vector<Point> points = array_to_points(mesh, verts);

            if (orient) {
                bool success = PMP::orient_polygon_soup(points, faces);
                if (!success) {
                    throw std::runtime_error("Polygon orientation failed");
                }
            }
            PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
            return mesh;
        }))
        .def("to_polygon_soup", [](const Mesh& mesh) {
            std::vector<Point> verts;
            std::vector<std::vector<size_t>> faces;
            PMP::polygon_mesh_to_polygon_soup(mesh, verts, faces);

            // convert points to arrays
            const size_t nv = mesh.number_of_vertices();
            py::array_t<double, py::array::c_style> verts_out({nv, ndims(mesh)});
            auto rv = verts_out.mutable_unchecked<2>();
            for (size_t i = 0; i < nv; i++) {
                Point p = verts[i];
                for (auto j = 0; j < ndims(mesh); j++) {
                    rv(i, j) = CGAL::to_double(p[j]);
                }
            }

            // ditto faces
            const size_t nf = mesh.number_of_faces();
            py::array_t<size_t, py::array::c_style> faces_out({nf, size_t(3)});
            auto rf = faces_out.mutable_unchecked<2>();
            for (size_t i = 0; i < nf; i++) {
                for (size_t j = 0; j < 3; j++) {
                    rf(i, j) = faces[i][j];
                }
            }
            return std::make_tuple(verts_out, faces_out);
        })
        .def("construct_points", &construct_points<Mesh, Point, F>)
        .def_property_readonly("is_valid", [](const Mesh& mesh) { return mesh.is_valid(false); })
        .def_property_readonly("n_vertices", [](const Mesh& mesh) { return mesh.number_of_vertices(); })
        .def_property_readonly("n_faces", [](const Mesh& mesh) { return mesh.number_of_faces(); })
        .def_property_readonly("n_edges", [](const Mesh& mesh) { return mesh.number_of_edges(); })

        .def_property_readonly("vertices", [](const Mesh& mesh) {
            std::vector<V> verts;
            verts.reserve(mesh.number_of_vertices());
            for (V v : mesh.vertices()) {
                verts.emplace_back(v);
            }
            return verts;
        })
        .def_property_readonly("faces", [](const Mesh& mesh) {
            std::vector<F> faces;
            faces.reserve(mesh.number_of_faces());
            for (F f : mesh.faces()) {
                faces.emplace_back(f);
            }
            return faces;
        })
        .def_property_readonly("edges", [](const Mesh& mesh) {
            std::vector<E> edges;
            edges.reserve(mesh.number_of_edges());
            for (E e : mesh.edges()) {
                edges.emplace_back(e);
            }
            return edges;
        })
        .def("edge_vertices", [](const Mesh& mesh, const std::vector<E>& edges) {
            std::map<V, size_t> vert_idxs;
            size_t vi = 0;
            for (V v : mesh.vertices()) {
                vert_idxs[v] = vi;
                vi++;
            }

            const size_t ne = edges.size();
            py::array_t<size_t, py::array::c_style> verts({ne, size_t(2)});
            auto r = verts.mutable_unchecked<2>();
            for (auto i = 0; i < ne; i++) {
                for (auto j = 0; j < 2; j++) {
                    r(i, j) = vert_idxs[mesh.vertex(edges[i], j)];
                }
            }

            return verts;
        })
    ;
}


void init_mesh(py::module &m) {
    py::module sub = m.def_submodule("surface_mesh");

    sub.def("triangulate_constrained_contour_pair", [](
        const py::array_t<double>& p_in, const py::array_t<double>& q_in,
        const std::vector<size_t>& pidx, const std::vector<size_t>& qidx, const size_t np0, const size_t nq0
    ) {
        Mesh3 mesh;  // todo gross
        std::vector<Point_3> p = array_to_points(mesh, p_in);
        std::vector<Point_3> q = array_to_points(mesh, q_in);
        return constrained_contour_pair_mesh(p, q, pidx, qidx, np0, nq0);
    });

    // Don't really understand how pybind11, typedefs, and templates interact here
    // But these serve as both Mesh3 and Mesh2 indices, so don't need to redefine them for Mesh2
    py::class_<V3>(sub, "Vertex");
    py::class_<F3>(sub, "Face");
    py::class_<E3>(sub, "Edge");
    py::class_<H3>(sub, "Halfedge");

    py::class_<AABB_Tree3>(sub, "AABB_Tree3");

    define_mesh<Mesh3, Point_3, V3, F3, E3, H3>(sub, "Mesh3")
        .def("aabb_tree", [](const Mesh3& mesh) {
            AABB_Tree3 tree;
            PMP::build_AABB_tree(mesh, tree);
            return tree;
        })
        .def("locate_points", [](const Mesh3& mesh, const AABB_Tree3& tree, const py::array_t<double>& points) {
            auto rp = points.unchecked<2>();
            size_t np = rp.shape(0);
            std::vector<F3> faces;
            faces.reserve(np);
            py::array_t<double, py::array::c_style> bary_coords({np, size_t(3)});

            auto rbc = bary_coords.mutable_unchecked<2>();

            for (size_t i = 0; i < np; i++) {
                Point_3 point = {rp(i, 0), rp(i, 1), rp(i, 2)};
                Loc3 loc = PMP::locate_with_AABB_tree(point, tree, mesh);
                faces.emplace_back(loc.first);

                for (size_t j = 0; j < 3; j++) {
                    rbc(i, j) = loc.second[j];
                }
            }

            return std::make_tuple(faces, bary_coords);
        })
        .def("remesh", [](Mesh3& mesh, const std::vector<F3>& faces, double target_edge_length, unsigned int n_iter) {
            PMP::isotropic_remeshing(faces, target_edge_length, mesh);
        })
        .def("remesh", [](Mesh3& mesh, const std::vector<F3>& faces, double target_edge_length, unsigned int n_iter, bool protect_constraints) {
            auto params = PMP::parameters::protect_constraints(protect_constraints);
            PMP::isotropic_remeshing(faces, target_edge_length, mesh, params);
        })
        .def("remesh", [](Mesh3& mesh, const std::vector<F3>& faces, double target_edge_length, unsigned int n_iter,
                          Mesh3::Property_map<V3, bool>& touched) {
            auto points = mesh.points();
            VertexPointMapWrapper point_map = VertexPointMapWrapper(points, touched);
            auto params = PMP::parameters::number_of_iterations(n_iter).vertex_point_map(point_map);
            PMP::isotropic_remeshing(faces, target_edge_length, mesh, params);
        })
        .def("expand_selection", [](Mesh3& mesh, const std::vector<V3>& selected) {
            std::set<V3> expanded;

            for (V3 v0 : selected) {
                expanded.insert(v0);

                for (V3 v1 : vertices_around_target(mesh.halfedge(v0), mesh)) {
                    expanded.insert(v1);
                }
            }

            return expanded;
        })
        .def("expand_selection", [](Mesh3& mesh, const std::vector<F3>& selected) {
            std::set<F3> expanded;

            for (F3 f0 : selected) {
                expanded.insert(f0);

                for (F3 f1 : faces_around_face(mesh.halfedge(f0), mesh)) {
                    expanded.insert(f1);
                }
            }

            return expanded;
        })
        .def("fair", [](Mesh3& mesh, const std::vector<V3>& verts, const py::kwargs& kwargs) {
            // A value controling the tangential continuity of the output surface patch.
            // The possible values are 0, 1 and 2, refering to the C0, C1 and C2 continuity.
            auto params = PMP::parameters::all_default();

            for (auto item : kwargs) {
                // auto key = item.first.cast<std::string>();
                auto key = py::cast<std::string>(item.first);

                if (key.compare("fairing_continuity")) {
                    //params.fairing_continuity(item.second.cast<unsigned int>());
                    params.fairing_continuity(py::cast<unsigned int>(item.second));
                }
            }

            // bool success = PMP::fair(mesh, verts, PMP::parameters::fairing_continuity(continuity));
            bool success = PMP::fair(mesh, verts, params);
            if (!success) {
                throw std::runtime_error("Fairing failed");
            }
        })
        .def("refine", [](Mesh3& mesh, const std::vector<F3>& faces, double density) {
            std::vector<V3> new_verts;
            std::vector<F3> new_faces;
            auto params = PMP::parameters::density_control_factor(density);
            PMP::refine(mesh, faces, std::back_inserter(new_faces), std::back_inserter(new_verts), params);
            return std::make_tuple(new_verts, new_faces);
        })
        .def("smooth_mesh", [](Mesh3& mesh, const std::vector<F3>& faces, unsigned int n_iter, bool use_safety_constraints) {
            auto params = PMP::parameters::number_of_iterations(n_iter).use_safety_constraints(use_safety_constraints);
            PMP::smooth_mesh(faces, mesh, params);
        })
        .def("smooth_shape", [](Mesh3& mesh, const std::vector<F3>& faces, const double time, unsigned int n_iter) {
            auto params = PMP::parameters::number_of_iterations(n_iter);
            PMP::smooth_shape(faces, mesh, time, params);
        })
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
            std::vector<Point_3> points;
            shortest_path.shortest_path_points_to_source_points(tgt_face, tgt_bc_, std::back_inserter(points));

            return points_to_array(points);
        })
        .def("face_normals", [](const Mesh3& mesh, const std::vector<F3> faces) {
            const size_t nf = faces.size();
            py::array_t<double, py::array::c_style> normals({nf, size_t(3)});
            auto r = normals.mutable_unchecked<2>();
            for (auto i = 0; i < nf; i++) {
                auto normal = PMP::compute_face_normal(faces[i], mesh);
                for (auto j = 0; j < 3; j++) {
                    r(i, j) = CGAL::to_double(normal[j]);
                }
            }
            return normals;
        })
    ;

   define_mesh<Mesh2, Point_2, V2, F2, E2, H2>(sub, "Mesh2")
//        .def("locate_with_aabb_tree", [](const Mesh2& mesh, const py::array_t<double>& points_in) {
//            using AABB_face_graph_primitive = typename CGAL::AABB_face_graph_triangle_primitive<Mesh2>;
//            using AABB_face_graph_traits = CGAL::AABB_traits<Kernel, AABB_face_graph_primitive>;
//            using AABBT = CGAL::AABB_tree<AABB_face_graph_traits>;
//            AABBT tree;
//            const Point2_to_Point3 vpm(mesh);
//            auto params = CGAL::parameters::vertex_point_map(vpm);
//            PMP::build_AABB_tree(mesh, tree, params);
//            const Point_3 pt = Point_3(0, 0, 0);
//            auto loc = PMP::locate_with_AABB_tree(pt, tree, mesh, params);
//            return true;
//
////            auto r = points_in.unchecked<2>();
////            size_t np = r.shape(0);
////            std::vector<Point_2> points;
////            points.reserve(np);
////
////            for (size_t i = 0; i < np; i++) {
////                points.emplace_back(Point_2(r(i, 0), r(i, 1)));
////            }
//
//            // return locate_points_with_aabb_tree<Mesh2, F2, Point_2, AABBT>(mesh, points, tree);
//            // auto loc = PMP::locate_with_AABB_tree(points[i], tree, mesh)
//        })
   ;
}
