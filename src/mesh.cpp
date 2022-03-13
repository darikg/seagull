#include "skgeom.hpp"
#include "funcs.hpp"

#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/refine.h>

#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Surface_mesh/IO/OFF.h>

typedef CGAL::Surface_mesh<Point_3>                                 Mesh;
typedef Mesh::Vertex_index                                          V;
typedef Mesh::Face_index                                            F;
typedef Mesh::Halfedge_index                                        H;
typedef Mesh::Edge_index                                            E;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename Key, typename Val>
Mesh::Property_map<Key, Val> add_property_map (Mesh& mesh, std::string name, const Val default_val) {
    Mesh::Property_map<Key, Val> pmap;
    bool created;
    std::tie(pmap, created) = mesh.add_property_map<Key, Val>(name, default_val);
    if (!created) {
        throw std::runtime_error("Property map already exists");
    }
    return pmap;
}

template <typename Key, typename Val>
Mesh::Property_map<Key, Val> get_property_map (Mesh& mesh, const std::string& name) {
    Mesh::Property_map<Key, Val> pmap;
    bool found;
    std::tie(pmap, found) = mesh.property_map<Key, Val>(name);
    if (!found) {
        throw std::runtime_error("Property map " + name + " not found");
    }
    return pmap;
}

template <typename Key, typename Val>
Val get_property_value(const Mesh::Property_map<Key, Val>& pmap, const Key& key) {
    return pmap[key];
}

template <typename Key, typename Val>
py::array_t<Val> get_property_values(const Mesh::Property_map<Key, Val>& pmap, const std::vector<Key>& keys) {
    size_t nk = keys.size();
    py::array_t<Val, py::array::c_style> vals({int(nk)});
    auto r = vals.mutable_unchecked<1>();

    for (size_t i = 0; i < nk; i++) {
        r(i) = pmap[keys[i]];
    }
    return vals;
}

template <typename Key, typename Val>
void set_property_value(Mesh::Property_map<Key, Val>& pmap, const Key& key, const Val val) {
    pmap[key] = val;
}

template <typename Key, typename Val>
void set_property_values(
        const Mesh::Property_map<Key, Val>& pmap, const std::vector<Key>& keys, const std::vector<Val>& vals)
{
    size_t nk = keys.size();
    size_t nv = vals.size();
    if (nk != nv) {
        throw std::runtime_error("Key and value array sizes do not match");
    }
    for (size_t i = 0; i < nk; i++) {
        pmap[keys[i]] = vals[i];
    }
}

struct VertexPointMapWrapper {
    // Used for tracking which verts get moved during remesh, etc
    using key_type = V;
    using value_type = Point_3;
    using reference = Point_3&;
    using category = boost::read_write_property_map_tag;

    Mesh::Property_map<V, Point_3>& points;
    Mesh::Property_map<V, bool>& touched;

    VertexPointMapWrapper(Mesh::Property_map<V, Point_3>& p, Mesh::Property_map<V, bool>& t) : points(p), touched(t) {}

    friend Point_3 get (const VertexPointMapWrapper& map, V v) {return map.points[v];}
    friend void put (const VertexPointMapWrapper& map, V v, const Point_3& point) {
        map.points[v] = point;
        map.touched[v] = true;
    }
};

struct CorefinementVisitor : public PMP::Corefinement::Default_visitor<Mesh> {
    // Used for tracking for refinement indices
    // CGAL's corefine only uses a visitor for the first mesh, so we need the references to both
    // here to tell which is which
    Mesh& mesh1;
    Mesh& mesh2;
    Mesh::Property_map<V, ssize_t>& vert_ids1;
    Mesh::Property_map<V, ssize_t>& vert_ids2;

    CorefinementVisitor(
        Mesh& m1, Mesh& m2, Mesh::Property_map<V, ssize_t>& v1, Mesh::Property_map<V, ssize_t>& v2
    ) : mesh1(m1), mesh2(m2), vert_ids1(v1), vert_ids2(v2) {}

    void new_vertex_added(size_t i_id, V v, const Mesh& mesh) {
        // Called when a new vertex is added in the mesh
        // (either an edge split or a vertex inserted in the interior of a face).
        // i_id is the intersection point id reported in new_node_added.
        // For each mesh, a vertex with a given id will be reported exactly once,
        // except if it is already an existing vertex.
        if (&mesh == &mesh1) {
            vert_ids1[v] = size_t(i_id);
        } else {
            vert_ids2[v] = size_t(i_id);
        }
    }
};


void init_mesh(py::module &m) {
    py::module sub = m.def_submodule("surface_mesh");

    py::class_<V>(sub, "Vertex");
    py::class_<F>(sub, "Face");
    py::class_<E>(sub, "Edge");
    py::class_<H>(sub, "Halfedge");

    py::class_<Mesh::Property_map<V, bool>>(sub, "VertexBoolProperty")
        .def("__getitem__", [](const Mesh::Property_map<V, bool> pmap, const V& vert) {
            return get_property_value(pmap, vert);
        })
        .def("__getitem__", [](const Mesh::Property_map<V, bool> pmap, const std::vector<V>& verts) {
            return get_property_values(pmap, verts);
        })
        .def("__setitem__", [](Mesh::Property_map<V, bool>& pmap, const V& vert, const bool val) {
            set_property_value(pmap, vert, val);
        })
        .def("__setitem__", [](Mesh::Property_map<V, bool>& pmap, const std::vector<V>& verts, const std::vector<bool>& vals) {
            set_property_values(pmap, verts, vals);
        })
    ;

    py::class_<Mesh::Property_map<V, ssize_t>>(sub, "VertexIntProperty")
        .def("__getitem__", [](const Mesh::Property_map<V, ssize_t> pmap, const V& vert) {
            return get_property_value(pmap, vert);
        })
        .def("__getitem__", [](const Mesh::Property_map<V, ssize_t> pmap, const std::vector<V>& verts) {
            return get_property_values(pmap, verts);
        })
        .def("__setitem__", [](Mesh::Property_map<V, ssize_t>& pmap, const V& vert, const ssize_t val) {
            set_property_value(pmap, vert, val);
        })
        .def("__setitem__", [](Mesh::Property_map<V, ssize_t>& pmap, const std::vector<V>& verts, const std::vector<ssize_t>& vals) {
            set_property_values(pmap, verts, vals);
        })
    ;

    py::class_<Mesh::Property_map<E, bool>>(sub, "EdgeBoolProperty")
        .def("__getitem__", [](const Mesh::Property_map<E, bool> pmap, const E& edge) {
            return get_property_value(pmap, edge);
        })
        .def("__getitem__", [](const Mesh::Property_map<E, bool> pmap, const std::vector<E>& edges) {
            return get_property_values(pmap, edges);
        })
        .def("__setitem__", [](Mesh::Property_map<E, bool>& pmap, const E& edge, const bool val) {
            set_property_value(pmap, edge, val);
        })
        .def("__setitem__", [](Mesh::Property_map<E, bool>& pmap, const std::vector<E>& edges, const std::vector<bool>& vals) {
            set_property_values(pmap, edges, vals);
        })
    ;

    py::class_<Mesh>(sub, "Mesh")
        .def(py::init<>())
        .def(py::init([](const std::string& file) {
            Mesh mesh;
            if(!CGAL::IO::read_polygon_mesh(file, mesh)) {
                throw std::runtime_error("Failed to load mesh");
            }
            return mesh;
        }))
        .def(py::init([](
                py::array_t<double> &verts,
                std::vector<std::vector<size_t>>& faces,
                const bool orient
            ) {
            // Convert verts to Point_3
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

            if (orient) {
                bool success = PMP::orient_polygon_soup(points, faces);
                if (!success) {
                    throw std::runtime_error("Polygon orientation failed");
                }
            }

            Mesh mesh;
            PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
            return mesh;
        }))
        .def("to_polygon_soup", [](const Mesh& mesh) {
            std::vector<Point_3> verts;
            std::vector<std::vector<size_t>> faces;
            PMP::polygon_mesh_to_polygon_soup(mesh, verts, faces);

            // convert points to arrays
            const size_t nv = mesh.number_of_vertices();
            py::array_t<double, py::array::c_style> verts_out({nv, size_t(3)});
            auto rv = verts_out.mutable_unchecked<2>();
            for (size_t i = 0; i < nv; i++) {
                Point_3 p = verts[i];
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

        .def("add_vertex_property", &add_property_map<V, bool>)
        .def("get_vertex_property", &get_property_map<V, bool>)
        .def("add_vertex_property", &add_property_map<V, ssize_t>)
        .def("get_vertex_property", &get_property_map<V, ssize_t>)
        .def("add_edge_property", &add_property_map<E, bool>)
        .def("get_edge_property", &get_property_map<E, bool>)

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
            for (size_t i = 0; i < ne; i++) {
                for (size_t j = 0; j < 2; j++) {
                    r(i, j) = vert_idxs[mesh.vertex(edges[i], j)];
                }
            }

            return verts;
        })

        .def("corefine", [](Mesh& mesh1, Mesh& mesh2){
            PMP::corefine(mesh1, mesh2);
        })
        .def("corefine", [](
                Mesh& mesh1, Mesh::Property_map<V, ssize_t>& vert_ids1, Mesh::Property_map<E, bool> ecm1,
                Mesh& mesh2, Mesh::Property_map<V, ssize_t>& vert_ids2, Mesh::Property_map<E, bool> ecm2) {

            CorefinementVisitor visitor(mesh1, mesh2, vert_ids1, vert_ids2);
            auto params1 = PMP::parameters::visitor(visitor).edge_is_constrained_map(ecm1);
            auto params2 = PMP::parameters::edge_is_constrained_map(ecm2);
            PMP::corefine(mesh1, mesh2, params1, params2);
        })
        .def("difference", [](Mesh& mesh1, Mesh& mesh2) {
            Mesh result;
            bool success = PMP::corefine_and_compute_difference(mesh1, mesh2, result);
            if (!success) {
                throw std::runtime_error("Boolean operation failed.");
            }
            return result;
        })
        .def("union", [](Mesh& mesh1, Mesh& mesh2) {
            Mesh result;
            bool success = PMP::corefine_and_compute_union(mesh1, mesh2, result);
            if (!success) {
                throw std::runtime_error("Boolean operation failed.");
            }
            return result;
        })
        .def("union", [](
                Mesh& mesh1, Mesh::Property_map<V, ssize_t>& vert_ids1, Mesh::Property_map<E, bool> ecm1,
                Mesh& mesh2, Mesh::Property_map<V, ssize_t>& vert_ids2, Mesh::Property_map<E, bool> ecm2) {

            CorefinementVisitor visitor(mesh1, mesh2, vert_ids1, vert_ids2);
            auto params1 = PMP::parameters::visitor(visitor).edge_is_constrained_map(ecm1);
            auto params2 = PMP::parameters::edge_is_constrained_map(ecm2);
            bool success = PMP::corefine_and_compute_union(mesh1, mesh2, mesh1);
            if (!success) {
                throw std::runtime_error("Boolean operation failed.");
            }
        })
        .def("intersection", [](Mesh& mesh1, Mesh& mesh2) {
            Mesh result;                        //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! copy paste mistake
            bool success = CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh1, mesh2, result);
            if (!success) {
                throw std::runtime_error("Boolean operation failed.");
            }
            return result;
        })
        .def("write_ply", [](Mesh& mesh, std::string file) {
            std::ofstream out(file, std::ios::binary);
            CGAL::IO::set_binary_mode(out);
            bool success = CGAL::IO::write_PLY(out, mesh, "", PMP::parameters::all_default());
            if (!success) {
                throw std::runtime_error("writing failed");
            }
        })
        .def("write_off", [](Mesh& mesh, std::string file) {
            std::ofstream out(file);
            bool success = CGAL::IO::write_OFF(out, mesh, PMP::parameters::all_default());
            if (!success) {
                throw std::runtime_error("writing failed");
            }
        })

        .def("remesh", [](Mesh& mesh, const std::vector<F>& faces, double target_edge_length, unsigned int n_iter) {
            PMP::isotropic_remeshing(faces, target_edge_length, mesh);
        })
        .def("remesh", [](Mesh& mesh, const std::vector<F>& faces, double target_edge_length, unsigned int n_iter,
                          Mesh::Property_map<V, bool>& touched) {
            auto points = mesh.points();
            VertexPointMapWrapper point_map = VertexPointMapWrapper(points, touched);
            auto params = PMP::parameters::number_of_iterations(n_iter).vertex_point_map(point_map);
            PMP::isotropic_remeshing(faces, target_edge_length, mesh, params);
        })
        .def("fair", [](Mesh& mesh, const std::vector<V>& verts, unsigned int continuity) {
            // A value controling the tangential continuity of the output surface patch.
            // The possible values are 0, 1 and 2, refering to the C0, C1 and C2 continuity.
            bool success = PMP::fair(mesh, verts, PMP::parameters::fairing_continuity(continuity));
            if (!success) {
                throw std::runtime_error("Fairing failed");
            }
        })
        .def("refine", [](Mesh& mesh, const std::vector<F>& faces, double density) {
            std::vector<V> new_verts;
            std::vector<F> new_faces;
            auto params = PMP::parameters::density_control_factor(density);
            PMP::refine(mesh, faces, std::back_inserter(new_faces), std::back_inserter(new_verts), params);
            return std::make_tuple(new_verts, new_faces);
        })
        .def("smooth_mesh", [](Mesh& mesh, const std::vector<F>& faces, unsigned int n_iter, bool use_safety_constraints) {
            auto params = PMP::parameters::number_of_iterations(n_iter).use_safety_constraints(use_safety_constraints);
            PMP::smooth_mesh(faces, mesh, params);
        })
    ;
}
