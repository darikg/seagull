#include "seagullmesh.hpp"


template <typename Mesh, typename Key, typename Val>
auto add_property_map(Mesh& mesh, std::string name, const Val default_val) {
    typename Mesh::Property_map<Key, Val> pmap;
    bool created;
    std::tie(pmap, created) = mesh.add_property_map<Key, Val>(name, default_val);
    if (!created) {
        throw std::runtime_error("Property map already exists");
    }
    return pmap;
}

template <typename Mesh, typename Key, typename Val>
auto define_property_map(py::module &m, std::string name) {
    // https://stackoverflow.com/a/47749076/7519203
    using PMap = typename Mesh::Property_map<Key, Val>;
    return py::class_<PMap>(m, name.c_str(), py::buffer_protocol(), py::dynamic_attr())
        .def("__getitem__", [](const PMap& pmap, const Key& key) {
            return pmap[key];
        })
        .def("__getitem__", [](const PMap& pmap, const std::vector<Key>& keys) {
            size_t nk = keys.size();
            py::array_t<Val, py::array::c_style> vals({int(nk)});
            auto r = vals.mutable_unchecked<1>();

            for (size_t i = 0; i < nk; i++) {
                r(i) = pmap[keys[i]];
            }
            return vals;
        })
        .def("__setitem__", [](PMap& pmap, const Key& key, const Val val) {
            pmap[key] = val;
        })
        .def("__setitem__", [](PMap& pmap, const std::vector<Key>& keys, const std::vector<Val>& vals) {
            size_t nk = keys.size();
            size_t nv = vals.size();
            if (nk != nv) {
                throw std::runtime_error("Key and value array sizes do not match");
            }
            for (size_t i = 0; i < nk; i++) {
                pmap[keys[i]] = vals[i];
            }
        })
    ;
}


template <typename Mesh, typename Key, typename Val>
void define_array_3_property_map(py::module &m, std::string name) {
    using PMap = typename Mesh::Property_map<Key, Val>;

    define_property_map<Mesh, Key, Val>(m, name)
        .def("get_array", [](const PMap& pmap, const std::vector<Key>& keys) {
            const size_t nk = keys.size();
            py::array_t<double, py::array::c_style> vals({nk, size_t(3)});
            auto r = vals.mutable_unchecked<2>();

            for (auto i = 0; i < nk; i++) {
                auto val = pmap[keys[i]];
                for (auto j = 0; j < 3; j++) {
                    r(i, j) = val[j];
                }
            }
            return vals;
        })
        .def("set_array", [](PMap& pmap, const std::vector<Key>& keys, const py::array_t<double>& vals) {
            const size_t nk = keys.size();
            auto r = vals.unchecked<2>();
            if (nk != r.shape(0)) {
                throw std::runtime_error("Key and value array sizes do not match");
            }
            if (3 != r.shape(1)) {
                throw std::runtime_error("Expected an array with 3 columns");
            }
            for (auto i = 0; i < nk; i++) {
                pmap[keys[i]] = Val(r(i, 0), r(i, 1), r(i, 2));
            }
        })
    ;
}

template <typename Mesh, typename Key, typename Val>
void define_array_2_property_map(py::module &m, std::string name) {
    using PMap = typename Mesh::Property_map<Key, Val>;

    define_property_map<Mesh, Key, Val>(m, name)
        .def("get_array", [](const PMap& pmap, const std::vector<Key>& keys) {
            const size_t nk = keys.size();
            py::array_t<double, py::array::c_style> vals({nk, size_t(2)});
            auto r = vals.mutable_unchecked<2>();

            for (auto i = 0; i < nk; i++) {
                auto val = pmap[keys[i]];
                for (auto j = 0; j < 2; j++) {
                    r(i, j) = val[j];
                }
            }
            return vals;
        })
        .def("set_array", [](PMap& pmap, const std::vector<Key>& keys, const py::array_t<double>& vals) {
            const size_t nk = keys.size();
            auto r = vals.unchecked<2>();
            if (nk != r.shape(0)) {
                throw std::runtime_error("Key and value array sizes do not match");
            }
            if (2 != r.shape(1)) {
                throw std::runtime_error("Expected an array with 2 columns");
            }
            for (auto i = 0; i < nk; i++) {
                pmap[keys[i]] = Val(r(i, 0), r(i, 1));
            }
        })
    ;
}

template<typename Mesh, typename V, typename F, typename E, typename H>
void define_mesh_properties(py::module &m, std::string name) {
    define_property_map        < Mesh, V, bool     >(m, "VertBoolPropertyMap"    + name);
    define_property_map        < Mesh, V, ssize_t  >(m, "VertIntPropertyMap"     + name);
    define_array_3_property_map< Mesh, V, Point3   >(m, "VertPoint3PropertyMap"  + name);
    define_array_3_property_map< Mesh, V, Vector3  >(m, "VertVector3PropertyMap" + name);
    define_array_2_property_map< Mesh, V, Point2   >(m, "VertPoint2PropertyMap"  + name);
    define_array_2_property_map< Mesh, V, Vector2  >(m, "VertVector2PropertyMap" + name);

    define_property_map<         Mesh, F, bool     >(m, "FaceBoolPropertyMap"    + name);
    define_property_map<         Mesh, F, ssize_t  >(m, "FaceIntPropertyMap"     + name);
    define_property_map<         Mesh, E, bool     >(m, "EdgeBoolPropertyMap"    + name);



    m.def("add_vertex_property", &add_property_map<Mesh, V, bool>)
     .def("add_vertex_property", &add_property_map<Mesh, V, ssize_t>)

     .def("add_vertex_property", &add_property_map<Mesh, V, Point3>)
     .def("add_vertex_property", &add_property_map<Mesh, V, Vector3>)
     .def("add_vertex_property", &add_property_map<Mesh, V, Point2>)
     .def("add_vertex_property", &add_property_map<Mesh, V, Vector2>)

     .def("add_face_property", &add_property_map<Mesh, F, bool>)
     .def("add_face_property", &add_property_map<Mesh, F, ssize_t>)

     .def("add_edge_property", &add_property_map<Mesh, E, bool>);
}

void init_properties(py::module &m) {
    py::module sub = m.def_submodule("properties");
    define_mesh_properties<Mesh3, V3, F3, E3, H3>(sub, "3");
    // define_mesh_properties<Mesh2, V2, F2, E2, H2>(sub, "2");
}