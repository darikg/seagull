from __future__ import annotations
from abc import ABC
from typing import Generic, TypeVar, Any, Optional, Dict, TYPE_CHECKING

from numpy import ndarray, zeros_like, array, sqrt, concatenate, repeat, ones

from seagullmesh._seagullmesh.mesh import (  # noqa
    Mesh3 as _Mesh3,
    polygon_soup_to_mesh3,
    Point2
    # Point3
)
from seagullmesh import _seagullmesh as sgm

if TYPE_CHECKING:
    import pyvista as pv  # noqa


class Mesh3:
    def __init__(self, mesh: Mesh3):
        self._mesh = mesh
        self.vertex_data = MeshData(mesh, sgm.properties.add_vertex_property, 'vertices')
        self.edge_data = MeshData(mesh, sgm.properties.add_edge_property, 'edges')

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

    def edge_vertices(self, edges):
        return self._mesh.edge_vertices(edges)

    def expand_selection(self, selection):
        return self._mesh.expand_selection(selection)

    def to_polygon_soup(self):
        return self._mesh.to_polygon_soup()

    def face_normals(self, faces):
        return self._mesh.face_normals(faces)

    @staticmethod
    def from_polygon_soup(verts: ndarray, faces: ndarray, orient=True) -> Mesh3:
        mesh = polygon_soup_to_mesh3(verts, faces, orient)
        return Mesh3(mesh)

    def to_pyvista(self) -> pv.PolyData:
        from pyvista import PolyData  # noqa
        verts, _faces = self._mesh.to_polygon_soup()
        faces = concatenate([3 * ones((_faces.shape[0], 1), dtype='int'), _faces.astype('int')], axis=1)
        mesh = PolyData(verts, faces=faces.reshape(-1))
        return mesh

    @staticmethod
    def from_pyvista(polydata: pv.PolyData, orient=True) -> Mesh3:
        assert polydata.is_all_triangles
        from pyvista._vtk import vtk_to_numpy  # noqa
        cells = vtk_to_numpy(polydata.GetPolys().GetConnectivityArray())
        faces = cells.reshape(-1, 3)
        return Mesh3.from_polygon_soup(polydata.points, faces, orient=orient)

    def corefine_tracked(self, other: Mesh3, vert_idx: str, edge_constrained: str):
        tracker, ecm1, ecm2 = _get_corefined_properties(self, other, vert_idx, edge_constrained)
        sgm.corefine.corefine(self._mesh, other._mesh, ecm1, ecm2, tracker)

    def union_tracked(self, other: Mesh3, vert_idx: str, edge_constrained: str):
        tracker, ecm1, ecm2 = _get_corefined_properties(self, other, vert_idx, edge_constrained)
        sgm.corefine.union(self._mesh, other._mesh, ecm1, ecm2, tracker)

    def remesh(
            self,
            faces,
            target_edge_length: float,
            n_iter: int,
            protect_constraints=False,
            touched_map: Optional[str] = None,
    ):
        if touched_map:
            touched = self.vertex_data.get_or_create_property(touched_map, default=False)
            sgm.meshing.remesh(self._mesh, faces, target_edge_length, n_iter, protect_constraints, touched)
        else:
            sgm.meshing.remesh(self._mesh, faces, target_edge_length, n_iter, protect_constraints)

    def fair(self, verts, continuity=0):
        sgm.meshing.fair(self._mesh, verts, continuity)

    def refine(self, faces, density=sqrt(3)):
        return sgm.meshing.refine(self._mesh, faces, density)

    def aabb_tree(self, points: Optional[str] = None):
        if points:
            return sgm.locate.aabb_tree(self._mesh, self.vertex_data[points])
        else:
            return self._mesh.aabb_tree()

    def locate_points(self, points: ndarray, aabb_tree=None):
        tree = aabb_tree or self.aabb_tree()
        return self._mesh.locate_points(tree, points)

    def construct_points(self, faces, bary_coords: ndarray) -> ndarray:
        return self._mesh.construct_points(faces, bary_coords)

    def lscm(self, uv_map: str):
        uv_map = self.vertex_data.get_or_create_property(uv_map, default=Point2(0, 0))
        sgm.parametrize.lscm(self._mesh, uv_map)

    def arap(self, uv_map: str):
        uv_map = self.vertex_data.get_or_create_property(uv_map, default=Point2(0, 0))
        sgm.parametrize.arap(self._mesh, uv_map)


def _get_corefined_properties(mesh1: Mesh3, mesh2: Mesh3, vert_idx: str, edge_constrained: str):
    vert_idx1 = mesh1.vertex_data.get_or_create_property(vert_idx, default=-1)
    vert_idx2 = mesh2.vertex_data.get_or_create_property(vert_idx, default=-1)
    tracker = sgm.CorefinementVertexTracker(mesh1.mesh, mesh2.mesh, vert_idx1, vert_idx2)
    ecm1 = mesh1.edge_data.get_or_create_property(edge_constrained, default=False)
    ecm2 = mesh2.edge_data.get_or_create_property(edge_constrained, default=False)
    return tracker, ecm1, ecm2


class MeshData:
    def __init__(self, mesh: Mesh3, add_fn, key_name: str):
        self._data = {}
        self._mesh = mesh
        self._add_fn = add_fn
        self._key_name = key_name

    def add_property(self, key: str, default: Any):
        pmap = self._add_fn(self._mesh, key, default)
        self._data[key] = pmap
        return pmap

    def get_or_create_property(self, key: str, default: Optional[Any] = None):
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
        keys = getattr(self._mesh, self._key_name)
        pmap[keys] = value
