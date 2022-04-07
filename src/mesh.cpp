#include "seagullmesh.hpp"
#include "util.hpp"

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

namespace PMP = CGAL::Polygon_mesh_processing;

template<typename T>
auto define_simple_type_3(py::module &m, std::string name) {
    py::class_<T>(m, name.c_str(), py::module_local())
        .def(py::init<double, double, double>())
    ;
}

template<typename T>
auto define_simple_type_2(py::module &m, std::string name) {
    py::class_<T>(m, name.c_str(), py::module_local())
        .def(py::init<double, double>())
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
                const bool orient
        ) {
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
        })
        .def("load_mesh_from_file", [](const std::string filename) {
            Mesh3 mesh;
            if(!CGAL::IO::read_polygon_mesh(filename, mesh)) {
                throw std::runtime_error("Failed to load mesh");
            }
            return mesh;
        })
    ;


    py::class_<V3>(sub, "Vertex");
    py::class_<F3>(sub, "Face");
    py::class_<E3>(sub, "Edge");
    py::class_<H3>(sub, "Halfedge");

    py::class_<Mesh3>(sub, "Mesh3")
        .def(py::init<>())
        
        .def_property_readonly("is_valid", [](const Mesh3& mesh) { return mesh.is_valid(false); })
        .def_property_readonly("n_vertices", [](const Mesh3& mesh) { return mesh.number_of_vertices(); })
        .def_property_readonly("n_faces", [](const Mesh3& mesh) { return mesh.number_of_faces(); })
        .def_property_readonly("n_edges", [](const Mesh3& mesh) { return mesh.number_of_edges(); })
        .def_property_readonly("n_halfedges", [](const Mesh3& mesh) { return mesh.number_of_halfedges(); })
        .def_property_readonly("points", [](const Mesh3& mesh) { return mesh.points(); })

        .def_property_readonly("vertices", [](const Mesh3& mesh) {
            std::vector<V3> verts;
            verts.reserve(mesh.number_of_vertices());
            for (V3 v : mesh.vertices()) {
                verts.emplace_back(v);
            }
            return verts;
        })
        .def_property_readonly("faces", [](const Mesh3& mesh) {
            std::vector<F3> faces;
            faces.reserve(mesh.number_of_faces());
            for (F3 f : mesh.faces()) {
                faces.emplace_back(f);
            }
            return faces;
        })
        .def_property_readonly("edges", [](const Mesh3& mesh) {
            std::vector<E3> edges;
            edges.reserve(mesh.number_of_edges());
            for (E3 e : mesh.edges()) {
                edges.emplace_back(e);
            }
            return edges;
        })
        .def_property_readonly("halfedges", [](const Mesh3& mesh) {
            std::vector<H3> halfedges;
            halfedges.reserve(mesh.number_of_halfedges());
            for (H3 h : mesh.halfedges()) {
                halfedges.emplace_back(h);
            }
            return halfedges;
        })

        .def("edge_vertices", [](const Mesh3& mesh, const std::vector<E3>& edges) {
            std::map<V3, size_t> vert_idxs;
            size_t vi = 0;
            for (V3 v : mesh.vertices()) {
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
        .def("to_polygon_soup", [](const Mesh3& mesh) {
            std::vector<Point3> verts;
            std::vector<std::vector<size_t>> faces;
            PMP::polygon_mesh_to_polygon_soup(mesh, verts, faces);

            // convert points to arrays
            // todo use util
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
        .def("estimate_geodesic_distances", [](const Mesh3& mesh, Mesh3::Property_map<V3, double>& distances, V3 source) {
            CGAL::Heat_method_3::estimate_geodesic_distances(mesh, distances, source);
        })
        .def("estimate_geodesic_distances", [](
                const Mesh3& mesh, Mesh3::Property_map<V3, double>& distances, const std::vector<V3>& sources) {
            CGAL::Heat_method_3::estimate_geodesic_distances(mesh, distances, sources);
        })

        .def("write_ply", [](Mesh3& mesh, std::string file) {
            std::ofstream out(file, std::ios::binary);
            CGAL::IO::set_binary_mode(out);
            bool success = CGAL::IO::write_PLY(out, mesh, "");
            if (!success) {
                throw std::runtime_error("writing failed");
            }
        })
        .def("write_off", [](Mesh3& mesh, std::string file) {
            std::ofstream out(file);
            bool success = CGAL::IO::write_OFF(out, mesh);
            if (!success) {
                throw std::runtime_error("writing failed");
            }
        })
    ;

   //define_mesh<Mesh2, Point2, V2, F2, E2, H2>(sub, "Mesh2");
}
