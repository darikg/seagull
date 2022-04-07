# Seagull Mesh

Convenience bindings to CGAL's [Surface Mesh](https://doc.cgal.org/latest/Surface_mesh/index.html)
 and [Polygon Mesh Processing](https://doc.cgal.org/latest/Polygon_mesh_processing/index.html) modules, plus some other 
assorted extras.

## Installation

```shell
conda create --name seagull -c conda-forge cgal-cpp eigen pybind11 pyvista  # (pyvista optional)
conda activate seagull
```

On linux, you'll also need

```shell
conda install -c conda-forge gxx_linux-64 libgcc
```

Finally, 

```shell
git clone git@github.com:darikg/seagull.git
pip install -vv ./seagull 
```

## Usage

### Overview

This package provides bindings to the `Surface_mesh<Point_3>` class using `Exact_predicates_inexact_constructions_kernel`. Because CGAL is a heavily templated library, only a few convenience bindings per overloaded method are available, more added as-needed. To simplify some of the python/C++ mappings, the bound classes and methods are further wrapped in a python layer.

### Mesh IO

3 constructors are available:
```python
from seagullmesh import Mesh3

mesh = Mesh3.from_polygon_soup(verts, faces)
mesh = Mesh3.from_file('mesh.ply')
mesh = Mesh3.from_pyvista(polydata)
```

with the corresponding output methods `mesh.to_polygon_soup`, `mesh.to_file`, `mesh.to_pyvista`.

### Mesh indices

CGAL indices `Surface_mesh::vertex_index`, `::edge_index`, `::face_index`, `::edge_index`, and `::halfedge_index` are exposed in the python properties `mesh.vertices`, `.faces`, `.edges`, and `.halfedges`, which are returned as numpy arrays of indices for convenience. These arrays are used for indexing property maps and specifying regions for further processing. (See below.)

### Mesh property maps

Property maps are stored in the python properties `mesh.vertex_data`, `.edge_data`, etc., which acts as python dicts and are indexed with arrays of the indices described above.

```python
mesh.vertex_data['foo'] = np.arange(mesh.n_vertices)  # Creates a new property
mesh.vertex_data['foo'][mesh.vertices[:3]] = [10, 11, 12]  # sub-indexing
foo_vals = mesh.vertex_data['foo'][mesh.vertices]
property_map = mesh.vertex_data.add_property('bar', default=1)  # Create a property map manually
```

Non-scalar valued properties are supported in the form of CGAL's `Point_2`, `Point_3`, `Vector_2`, `Vector_3`objects. These property maps require some special handling for transforming to/from numpy arrays in the form of the`set_array`/`get_array` instead of direct indexing.

```python
from seagullmesh import Point2
mesh.vertex_data.add_property('uv_map', default=Point2(0, 0))
mesh.vertex_data['uv_map'].set_array(mesh.vertices, np.random.uniform(-1, 1, (mesh.n_vertices, 2)))
first_2_points = mesh.vertex_data['uv_map'].get_array(mesh.vertices[:2])
```

## Mesh processing

Currently implemented are:

From [PMP Meshing](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html):
  - `mesh.remesh(faces, target_edge_length, n_iterations)`
  - `mesh.fair(vertices, fairing_continuity)`
  - `mesh.refine(faces, density)`
  - `mesh.smooth_mesh(mesh, faces, n_iterations)`
  - `mesh.smooth_shape(mesh, faces, time)`

From [PMP Corefinement and Boolean Operations](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html)
  - `mesh.corefine(other)`
  - `mesh.union(other)`
  - `mesh.difference(other)`
  - `mesh.intersection(other)`

From the [PMP Location Functions](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__locate__grp.html)
  - `mesh.aabb_tree()`
  - `mesh.locate_points(points, aabb_tree)`
  - `mesh.construct_points(faces, bary_coords)`

From [Triangulated Surface Mesh Shortest Paths
](https://doc.cgal.org/latest/Surface_mesh_shortest_path/group__PkgSurfaceMeshShortestPathRef.html)
  - `mesh.shortest_path(src_face, src_bary_coords, tgt_face, tgt_bary_coords)`

From [The Heat Method](https://doc.cgal.org/latest/Heat_method_3/classCGAL_1_1Heat__method__3_1_1Surface__mesh__geodesic__distances__3.html)
  - `mesh.estimate_geodesic_distances(source_vertex_or_vertices)`

From [Planar Parameterization of Triangulated Surface Meshes](https://doc.cgal.org/latest/Surface_mesh_parameterization/group__PkgSurfaceMeshParameterizationRef.html)
  - `mesh.lscm(uv_map)`
  - `mesh.arap(uv_map)`

## Acknowledgements

The basis of the pybind11-cgal infrastructure is inspired very heavily by [Scikit-geometry](https://github.com/scikit-geometry/scikit-geometry).

## See also
  - [Scikit-geometry](https://github.com/scikit-geometry/scikit-geometry)
  - [pygalmesh](https://github.com/meshpro/pygalmesh)
  - [potpourri3d](https://github.com/nmwsharp/potpourri3d)
  - [pyvista](https://github.com/pyvista/pyvista)

## License

This software is licensed under the LGPL-3 license. See the [LICENSE](LICENSE) file for details.
