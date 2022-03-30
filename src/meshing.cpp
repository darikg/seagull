#include "skgeom.hpp"

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/refine.h>

namespace PMP = CGAL::Polygon_mesh_processing;

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

void init_meshing(py::module &m) {
    py::module sub = m.def_submodule("meshing");
    sub.def("remesh", [](Mesh3& mesh, const std::vector<F3>& faces, double target_edge_length, unsigned int n_iter) {
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
    ;
}
