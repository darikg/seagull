#include "skgeom.hpp"
#include <CGAL/Polygon_mesh_processing/corefinement.h>

namespace PMP = CGAL::Polygon_mesh_processing;

struct CorefinementVisitor : public PMP::Corefinement::Default_visitor<Mesh3> {
    // Used for tracking for refinement indices
    // CGAL's corefine only uses a visitor for the first mesh, so we need the references to both
    // here to tell which is which
    Mesh3& mesh1;
    Mesh3& mesh2;
    Mesh3::Property_map<V3, ssize_t>& vert_ids1;
    Mesh3::Property_map<V3, ssize_t>& vert_ids2;

    CorefinementVisitor(
        Mesh3& m1, Mesh3& m2, Mesh3::Property_map<V3, ssize_t>& v1, Mesh3::Property_map<V3, ssize_t>& v2
    ) : mesh1(m1), mesh2(m2), vert_ids1(v1), vert_ids2(v2) {}

    void new_vertex_added(size_t i_id, V3 v, const Mesh3& mesh) {
        // Called when a new vertex is added in the mesh
        // (either an edge split or a vertex inserted in the interior of a face).
        // i_id is the intersection point id reported in new_node_added.
        // For each mesh, a vertex with a given id will be reported exactly once,
        // except if it is already an existing vertex.
        if (&mesh == &mesh1) {
            vert_ids1[v] = ssize_t(i_id);
        } else {
            vert_ids2[v] = ssize_t(i_id);
        }
    }
};

template <typename Mesh, typename E>
auto named_parameters(const py::dict& kwargs) {
    using ECM = typename Mesh::Property_map<E, bool>;
    auto params = PMP::parameters::all_default();

    for (auto item : kwargs) {
        auto key = py::cast<std::string>(item.first);

        if (key.compare("edge_is_constrained_map")) {
            params.edge_is_constrained_map(py::cast<ECM>(item.second));
        }
    }
    return params;
}

void init_corefine(py::module &m) {
    py::module sub = m.def_submodule("corefine");
    sub.def("corefine", [](Mesh3& mesh1, Mesh3& mesh2){
        PMP::corefine(mesh1, mesh2);
    })
    .def("corefine", [](Mesh3& mesh1, Mesh3& mesh2, const py::dict& kwargs1){
        auto np1 = named_parameters<Mesh3, E3>(kwargs1);
        PMP::corefine(mesh1, mesh2, np1);
    })
    // .def("corefine", [](
    //         Mesh3& mesh1, Mesh3::Property_map<V3, ssize_t>& vert_ids1, Mesh3::Property_map<E3, bool> ecm1,
    //         Mesh3& mesh2, Mesh3::Property_map<V3, ssize_t>& vert_ids2, Mesh3::Property_map<E3, bool> ecm2) {

    //     CorefinementVisitor visitor(mesh1, mesh2, vert_ids1, vert_ids2);
    //     auto params1 = PMP::parameters::visitor(visitor).edge_is_constrained_map(ecm1);
    //     auto params2 = PMP::parameters::edge_is_constrained_map(ecm2);
    //     PMP::corefine(mesh1, mesh2, params1, params2);
    // })
    // .def("difference", [](Mesh3& mesh1, Mesh3& mesh2) {
    //     Mesh3 result;
    //     bool success = PMP::corefine_and_compute_difference(mesh1, mesh2, result);
    //     if (!success) {
    //         throw std::runtime_error("Boolean operation failed.");
    //     }
    //     return result;
    // })
    // .def("union", [](Mesh3& mesh1, Mesh3& mesh2) {
    //     Mesh3 result;
    //     bool success = PMP::corefine_and_compute_union(mesh1, mesh2, result);
    //     if (!success) {
    //         throw std::runtime_error("Boolean operation failed.");
    //     }
    //     return result;
    // })
    // .def("union", [](
    //         Mesh3& mesh1, Mesh3::Property_map<V3, ssize_t>& vert_ids1, Mesh3::Property_map<E3, bool> ecm1,
    //         Mesh3& mesh2, Mesh3::Property_map<V3, ssize_t>& vert_ids2, Mesh3::Property_map<E3, bool> ecm2) {

    //     CorefinementVisitor visitor(mesh1, mesh2, vert_ids1, vert_ids2);
    //     auto params1 = PMP::parameters::visitor(visitor).edge_is_constrained_map(ecm1);
    //     auto params2 = PMP::parameters::edge_is_constrained_map(ecm2);
    //     bool success = PMP::corefine_and_compute_union(mesh1, mesh2, mesh1, params1, params2);
    //     if (!success) {
    //         throw std::runtime_error("Boolean operation failed.");
    //     }
    // })
    // .def("intersection", [](Mesh3& mesh1, Mesh3& mesh2) {
    //     Mesh3 result;
    //     bool success = CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh1, mesh2, result);
    //     if (!success) {
    //         throw std::runtime_error("Boolean operation failed.");
    //     }
    //     return result;
    // })
    ;
}