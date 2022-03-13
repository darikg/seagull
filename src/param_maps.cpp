// https://stackoverflow.com/a/47749076/7519203

template<typename Key, typename Val>
void define_property_map(py::module& m, std::string& kvtype) {
    using Class = Mesh::Property_Map<Key, Val>;
    std::string pyclass_name = kvtype + std::string("PropertyMap")
    py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
    .def(py::init<>())
    .def(py::init<Class::xy_t, Class::xy_t, T>())
    .def("size",      &Class::size)
    .def("width",     &Class::width)
    .def("height",    &Class::height);

    py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
        .def("__getitem__", [](const Class& pmap, const Key& key) {
            return get_property_value(pmap, key);
        })
        .def("__getitem__", [](const Class& pmap, const std::vector<Key>& keys) {
            return get_property_values(pmap, keys);
        })
        .def("__setitem__", [](Class& pmap, const Key& key, const Val val) {
            set_property_value(pmap, key, val);
        })
        .def("__setitem__", [](Class& pmap, const std::vector<Key>& keys, const std::vector<Val>& vals) {
            set_property_values(pmap, keys, vals);
        })
    ;
}