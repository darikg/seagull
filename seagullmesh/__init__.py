from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING, Union, Sequence, Protocol, TypeVar, overload, Tuple

from numpy import ndarray, zeros_like, array, sqrt, concatenate, ones
from seagullmesh._seagullmesh.mesh import (  # noqa
    Mesh3 as _Mesh3,
    polygon_soup_to_mesh3,
    load_mesh_from_file,
    Point2, Point3, Vector2, Vector3,
    Vertex, Face, Edge, Halfedge,
)

from seagullmesh import _seagullmesh as sgm
from ._version import version_info, __version__  # noqa

Vertices = Sequence[Vertex]
Faces = Sequence[Face]
Edges = Sequence[Edge]
Halfedges = Sequence[Halfedge]


if TYPE_CHECKING:
    try:
        import pyvista as pv  # noqa
    except ImportError:
        pass

A = ndarray


class Mesh3:
    def __init__(self, mesh: _Mesh3):
        self._mesh = mesh
        self.vertex_data = MeshData(mesh, sgm.properties.add_vertex_property, 'vertices')
        self.face_data = MeshData(mesh, sgm.properties.add_face_property, 'faces')
        self.edge_data = MeshData(mesh, sgm.properties.add_edge_property, 'edges')
        self.halfedge_data = MeshData(mesh, sgm.properties.add_halfedge_property, 'halfedges')

    mesh = property(lambda self: self._mesh)

    vertices = property(lambda self: array(self._mesh.vertices))
    faces = property(lambda self: array(self._mesh.faces))
    edges = property(lambda self: array(self._mesh.edges))
    halfedges = property(lambda self: array(self._mesh.halfedges))

    n_vertices = property(lambda self: array(self._mesh.n_vertices))
    n_faces = property(lambda self: array(self._mesh.n_faces))
    n_edges = property(lambda self: array(self._mesh.n_edges))
    n_halfedges = property(lambda self: array(self._mesh.n_halfedges))

    is_valid = property(lambda self: self._mesh.is_valid)
    points = property(lambda self: self._mesh.points)

    def edge_vertices(self, edges: Edges) -> A:
        """Returns a len(edges) * 2 array of integer vertex indices"""
        return self._mesh.edge_vertices(edges)

    def expand_selection(self, selection: Sequence[Key]) -> Sequence[Key]:
        """Given a list of vertices or faces, returns a sequence containing the original and adjacent elements"""
        return self._mesh.expand_selection(selection)

    def to_polygon_soup(self) -> Tuple[A, A]:
        """Returns vertices (nv * 3) and faces (nf * 3) array"""
        return self._mesh.to_polygon_soup()

    def face_normals(self, faces: Faces) -> A:
        """Returns a (len(faces) * 3) array of face normal vectors"""
        return self._mesh.face_normals(faces)

    @staticmethod
    def from_polygon_soup(verts: A, faces: A, orient=True) -> Mesh3:
        """Constructs a surface mesh from vertices (nv * 3) and faces (nf * 3) arrays

        If `orient` is True (default), the faces are reindexed to represent a consistent manifold surface.
        """
        mesh = polygon_soup_to_mesh3(verts, faces, orient)
        return Mesh3(mesh)

    @staticmethod
    def from_file(filename: str) -> Mesh3:
        mesh = load_mesh_from_file(filename)
        return Mesh3(mesh)

    def to_file(self, filename: str):
        ext = Path(filename).suffix
        if ext == '.ply':
            self._mesh.write_ply(filename)
        elif ext == '.off':
            self._mesh.write_off(filename)
        else:
            raise ValueError(f"Unsupported format '{ext}'")

    @staticmethod
    def from_pyvista(polydata: pv.PolyData, orient=True) -> Mesh3:
        """Converts a `pyvista.PolyData` object to a surface mesh.

        Currently, all point/cell data is ignored.
        """
        assert polydata.is_all_triangles
        from pyvista._vtk import vtk_to_numpy  # noqa
        cells = vtk_to_numpy(polydata.GetPolys().GetConnectivityArray())
        faces = cells.reshape(-1, 3)
        return Mesh3.from_polygon_soup(polydata.points, faces, orient=orient)

    def to_pyvista(self) -> pv.PolyData:
        """Returns the mesh as a `pyvista.PolyData` object.

        Currently, all vertex/face/edge/halfedge data is ignored.
        """
        from pyvista import PolyData  # noqa
        verts, _faces = self._mesh.to_polygon_soup()
        faces = concatenate([3 * ones((_faces.shape[0], 1), dtype='int'), _faces.astype('int')], axis=1)
        mesh = PolyData(verts, faces=faces.reshape(-1))
        return mesh

    def corefine(self, other: Mesh3) -> None:
        """Corefines the two meshes in place"""
        sgm.corefine.corefine(self._mesh, other._mesh)

    def union(self, other: Mesh3, inplace=False) -> Mesh3:
        """Corefines the two meshes and returns their boolean union"""
        out = self if inplace else Mesh3(_Mesh3())
        sgm.corefine.union(self._mesh, other._mesh, out._mesh)
        return out

    def difference(self, other: Mesh3, inplace=False) -> Mesh3:
        """Corefines the two meshes and returns their boolean difference"""
        out = self if inplace else Mesh3(_Mesh3())
        sgm.corefine.difference(self._mesh, other._mesh, out._mesh)
        return out

    def intersection(self, other: Mesh3, inplace=False) -> Mesh3:
        """Corefines the two meshes and returns their boolean intersection"""
        out = self if inplace else Mesh3(_Mesh3())
        sgm.corefine.intersection(self._mesh, other._mesh, out._mesh)
        return out

    def corefine_tracked(
            self,
            other: Mesh3,
            vert_idx: Union[str, PropertyMap[Vertex, int]],
            edge_constrained: Union[str, PropertyMap[Edge, bool]],
    ) -> None:
        tracker, ecm1, ecm2 = _get_corefined_properties(self, other, vert_idx, edge_constrained)
        sgm.corefine.corefine(self._mesh, other._mesh, ecm1, ecm2, tracker)

    def union_tracked(
            self,
            other: Mesh3,
            vert_idx: Union[str, PropertyMap[Vertex, int]],
            edge_constrained: Union[str, PropertyMap[Edge, bool]],
    ) -> None:
        tracker, ecm1, ecm2 = _get_corefined_properties(self, other, vert_idx, edge_constrained)
        sgm.corefine.union(self._mesh, other._mesh, ecm1, ecm2, tracker)

    def remesh(
            self,
            faces: Faces,
            target_edge_length: float,
            n_iter: int,
            protect_constraints=False,
            touched_map: Optional[Union[str, PropertyMap[Vertex, bool]]] = None,
    ) -> None:
        """Perform isotropic remeshing on the specified faces.

        If an optional `touched_map` mapping vertices to bools is specified, vertices that were created or moved during
        the remeshing are flagged as True.
        """
        if touched_map:
            touched = self.vertex_data.get_or_create_property(touched_map, default=False)
            sgm.meshing.remesh(self._mesh, faces, target_edge_length, n_iter, protect_constraints, touched)
        else:
            sgm.meshing.remesh(self._mesh, faces, target_edge_length, n_iter, protect_constraints)

    def fair(self, verts: Vertices, continuity=0) -> None:
        """Fair the specified mesh vertices"""
        sgm.meshing.fair(self._mesh, verts, continuity)

    def refine(self, faces: Faces, density=sqrt(3)) -> Tuple[Vertices, Faces]:
        """Refine the specified mesh faces

        The number of faces is increased by a factor of `density`.
        Returns indices to the newly created vertices and faces.
        """
        return sgm.meshing.refine(self._mesh, faces, density)

    def smooth_mesh(self, faces: Faces, n_iter: int, use_safety_constraints=False) -> None:
        """Smooth the specified faces"""
        sgm.meshing.smooth_mesh(self._mesh, faces, n_iter, use_safety_constraints)

    def smooth_shape(self, faces: Faces, time: float, n_iter: int) -> None:
        """Smooth the mesh shape by mean curvature flow

        A larger time step results in faster convergence but details may be distorted to a larger extent compared to
         more iterations with a smaller step. Typical values scale in the interval (1e-6, 1]
        """
        sgm.meshing.smooth_shape(self._mesh, faces, time, n_iter)

    def aabb_tree(self, vert_points: Optional[Union[str, VertexPointMap]] = None):
        """Construct an axis-aligned bounding box tree for accelerated point location by `Mesh3.locate_points

        By default, the AABB tree is constructed for the default mesh vertex locations, but also accepts a vertex
        property map storing Point2 or Point3 locations.
        """
        if vert_points:
            vert_points = self.vertex_data[vert_points] if isinstance(vert_points, str) else vert_points
            return sgm.locate.aabb_tree(self._mesh, vert_points)
        else:
            return sgm.locate.aabb_tree(self._mesh)

    def locate_points(self, points: ndarray, aabb_tree=None) -> Tuple[Faces, A]:
        """Given an array of points, locate the nearest corresponding points on the mesh

        `aabb_tree` is an optional axis-aligned bounding box from Mesh3.aabb_tree. If the tree was constructed with a
        Point2 vertex property map, `points` is of shape (np, 2), otherwise (np, 3).

        Returns a list of face indices of length np, and an array (np, 3) of barycentric coordinates within those faces.
        """
        tree = aabb_tree or self.aabb_tree()
        return sgm.locate.locate_points(self._mesh, tree, points)

    def construct_points(
            self,
            faces: Faces,
            bary_coords: A,
            vert_points: Optional[Union[str, VertexPointMap]] = None,
    ) -> A:
        """Construct a set of points from face barycentric coordinates

        `bary_coords` must be of shape (len(faces), 3)

        By default, points are constructed using the default mesh vertex points. An optional vertex point map of value
        Point2 or Point3 can also be supplied. The returned array if of shape (len(faces), 2 or 3) as appropriate.
        """
        if vert_points:
            vert_points = self.vertex_data[vert_points] if isinstance(vert_points, str) else vert_points
            return sgm.locate.construct_points(self._mesh, faces, bary_coords, vert_points)
        else:
            return sgm.locate.construct_points(self._mesh, faces, bary_coords)

    def shortest_path(
            self,
            src_face: Face,
            src_bc: A,
            tgt_face: Face,
            tgt_bc: A,
    ):
        """Constructs the shortest path between the source and target locations

        locations are specified as a face and barycentric coordinates
        """
        return sgm.locate.shortest_path(self._mesh, src_face, src_bc, tgt_face, tgt_bc)

    def lscm(self, uv_map: Union[str, UvMap]):
        """Performs least-squares conformal mapping"""
        uv_map = self.vertex_data.get_or_create_property(uv_map, default=Point2(0, 0))
        sgm.parametrize.lscm(self._mesh, uv_map)

    def arap(self, uv_map: Union[str, UvMap]):
        """Performs as-rigid-as-possible parameterization"""
        uv_map = self.vertex_data.get_or_create_property(uv_map, default=Point2(0, 0))
        sgm.parametrize.arap(self._mesh, uv_map)

    def estimate_geodesic_distances(
            self,
            src: Union[Vertex, Vertices],
            distance_prop: Union[str, PropertyMap[Vertex, float]],
    ):
        """Estimates the geodesic distance from the source vertex/vertices to all vertices in the mesh

        Estimated distances are stored in the supplied vertex property map.
        """
        distances = self.vertex_data.get_or_create_property(distance_prop, default=0.0)
        self._mesh.estimate_geodesic_distances(distances, src)


def _get_corefined_properties(mesh1: Mesh3, mesh2: Mesh3, vert_idx: str, edge_constrained: str):
    vert_idx1 = mesh1.vertex_data.get_or_create_property(vert_idx, default=-1)
    vert_idx2 = mesh2.vertex_data.get_or_create_property(vert_idx, default=-1)
    tracker = sgm.corefine.CorefinementVertexTracker(mesh1.mesh, mesh2.mesh, vert_idx1, vert_idx2)
    ecm1 = mesh1.edge_data.get_or_create_property(edge_constrained, default=False)
    ecm2 = mesh2.edge_data.get_or_create_property(edge_constrained, default=False)
    return tracker, ecm1, ecm2


Key = TypeVar('Key', Vertex, Face, Edge, Halfedge)
Val = TypeVar('Val', int, bool, Point2, Point3, Vector2, Vector3)


class PropertyMap(Protocol[Key, Val]):
    @overload
    def __getitem__(self, item: Key) -> Val: ...

    @overload
    def __getitem__(self, item: Sequence[Key]) -> Sequence[Val]: ...

    def __getitem__(self, item): ...

    @overload
    def __setitem__(self, key: Key, value: Val): ...

    @overload
    def __setitem__(self, key: Sequence[Key], value: Sequence[Val]): ...

    def __setitem__(self, key, item): ...


VertexPointMap = PropertyMap[Vertex, Union[Point2, Point3]]
UvMap = PropertyMap[Vertex, Point2]


class MeshData:
    def __init__(self, mesh: Mesh3, add_fn, key_name: str):
        self._data = {}
        self._mesh = mesh
        self._add_fn = add_fn
        self._key_name = key_name

    @property
    def mesh_keys(self):
        return getattr(self._mesh, self._key_name)

    @property
    def n_mesh_keys(self) -> int:
        return getattr(self._mesh, f'n_{self._key_name}')

    def add_property(self, key: str, default: Any):
        pmap = self._add_fn(self._mesh, key, default)
        self._data[key] = pmap
        return pmap

    def get_or_create_property(self, key: Union[str, PropertyMap], default: Optional[Any] = None):
        if not isinstance(key, str):
            return key  # assume it's a property map

        if key in self._data:
            return self._data[key]
        else:
            pmap = self._add_fn(self._mesh, key, default)
            self._data[key] = pmap
            return pmap

    def __getitem__(self, item: str):
        return self._data[item]

    def __setitem__(self, key: str, value: ndarray):
        default = zeros_like(value, shape=()).item()
        pmap = self.get_or_create_property(key, default)
        pmap[self.mesh_keys] = value
