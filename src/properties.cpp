#include "skgeom.hpp"


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
void define_property_map(py::module &m, std::string name) {
    // https://stackoverflow.com/a/47749076/7519203
    using PMap = typename Mesh::Property_map<Key, Val>;
    py::class_<PMap>(m, name.c_str(), py::buffer_protocol(), py::dynamic_attr())
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

template<typename Mesh, typename V, typename F, typename E, typename H>
void define_mesh_properties(py::module &m, std::string name) {
    define_property_map<Mesh, V, bool>(m, "VertBoolPropertyMap" + name);
    define_property_map<Mesh, V, ssize_t>(m, "VertIntPropertyMap" + name);
    define_property_map<Mesh, F, bool>(m, "FaceBoolPropertyMap" + name);
    define_property_map<Mesh, F, ssize_t>(m, "FaceIntPropertyMap" + name);
    define_property_map<Mesh, E, bool>(m, "EdgeBoolPropertyMap" + name);

    m.def("add_vertex_property", &add_property_map<Mesh, V, bool>)
     .def("add_vertex_property", &add_property_map<Mesh, V, ssize_t>)
     .def("add_face_property", &add_property_map<Mesh, F, bool>)
     .def("add_face_property", &add_property_map<Mesh, F, ssize_t>)
     .def("add_edge_property", &add_property_map<Mesh, E, bool>);
}

void init_properties(py::module &m) {
    py::module sub = m.def_submodule("properties");
    define_mesh_properties<Mesh3, V3, F3, E3, H3>(sub, "3");
    define_mesh_properties<Mesh2, V2, F2, E2, H2>(sub, "2");
}