#include "seagullmesh.hpp"

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Surface_mesh/IO/OFF.h>

namespace PMP = CGAL::Polygon_mesh_processing;


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

template<typename T>
auto define_simple_type_3(py::module &m, std::string name) {
    py::class_<T>(m, name.c_str())
        .def(py::init<>([](double x, double y, double z) {
            return T(x, y, z);
        }));
}

template<typename T>
auto define_simple_type_2(py::module &m, std::string name) {
    py::class_<T>(m, name.c_str())
        .def(py::init<>([](double x, double y) {
            return T(x, y);
        }));
}


template<typename Mesh, typename Point, typename V, typename F, typename E, typename H>
auto define_mesh(py::module &m, std::string name) {
    return py::class_<Mesh>(m, name.c_str())
        .def(py::init<>())
        
        .def_property_readonly("is_valid", [](const Mesh& mesh) { return mesh.is_valid(false); })
        .def_property_readonly("n_vertices", [](const Mesh& mesh) { return mesh.number_of_vertices(); })
        .def_property_readonly("n_faces", [](const Mesh& mesh) { return mesh.number_of_faces(); })
        .def_property_readonly("n_edges", [](const Mesh& mesh) { return mesh.number_of_edges(); })
        .def_property_readonly("n_halfedges", [](const Mesh& mesh) { return mesh.number_of_halfedges(); })
        .def_property_readonly("points", [](const Mesh& mesh) { return mesh.points(); })

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
        .def_property_readonly("halfedges", [](const Mesh& mesh) {
            std::vector<H> halfedges;
            halfedges.reserve(mesh.number_of_halfedges());
            for (H h : mesh.halfedges()) {
                halfedges.emplace_back(h);
            }
            return halfedges;
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
        .def("expand_selection", [](Mesh& mesh, const std::vector<V>& selected) {
            std::set<V> expanded;
            for (V v0 : selected) {
                expanded.insert(v0);
                for (V v1 : vertices_around_target(mesh.halfedge(v0), mesh)) {
                    expanded.insert(v1);
                }
            }
            return expanded;
        })
        .def("expand_selection", [](Mesh& mesh, const std::vector<F>& selected) {
            std::set<F> expanded;
            for (F f0 : selected) {
                expanded.insert(f0);
                for (F f1 : faces_around_face(mesh.halfedge(f0), mesh)) {
                    expanded.insert(f1);
                }
            }
            return expanded;
        })
    ;
}


void init_mesh(py::module &m) {
    py::module sub = m.def_submodule("mesh");

    define_simple_type_2<Point2>(sub, "Point2");
    define_simple_type_3<Point3>(sub, "Point3");
    define_simple_type_2<Vector2>(sub, "Vector2");
    define_simple_type_3<Vector3>(sub, "Vector3");

    sub.def("polygon_soup_to_mesh3", [](
            py::array_t<double> &points, 
            std::vector<std::vector<size_t>>& faces, 
            const bool orient) {

        Mesh3 mesh;
        std::vector<Point3> vertices = array_to_points_3(points);

        if (orient) {
            bool success = PMP::orient_polygon_soup(vertices, faces);
            if (!success) {
                throw std::runtime_error("Polygon orientation failed");
            }
        }
        PMP::polygon_soup_to_polygon_mesh(vertices, faces, mesh);
        return mesh;
    });

    // Don't really understand how pybind11, typedefs, and templates interact here
    // But these serve as both Mesh3 and Mesh2 indices, so don't need to redefine them for Mesh2
    py::class_<V3>(sub, "Vertex");
    py::class_<F3>(sub, "Face");
    py::class_<E3>(sub, "Edge");
    py::class_<H3>(sub, "Halfedge");

    define_mesh<Mesh3, Point3, V3, F3, E3, H3>(sub, "Mesh3")
        .def("to_polygon_soup", [](const Mesh3& mesh) {
            std::vector<Point3> verts;
            std::vector<std::vector<size_t>> faces;
            PMP::polygon_mesh_to_polygon_soup(mesh, verts, faces);

            // convert points to arrays
            const size_t nv = mesh.number_of_vertices();
            py::array_t<double, py::array::c_style> verts_out({nv, size_t(3)});
            auto rv = verts_out.mutable_unchecked<2>();
            for (size_t i = 0; i < nv; i++) {
                Point3 p = verts[i];
                for (auto j = 0; j < 3; j++) {
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

   //define_mesh<Mesh2, Point2, V2, F2, E2, H2>(sub, "Mesh2");
}
