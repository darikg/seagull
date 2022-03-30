from __future__ import annotations
from abc import ABC
from re import A
from typing import Generic, TypeVar, Any

from numpy import ndarray, zeros_like, array
from skgeom._skgeom.mesh import Mesh2, Mesh3
from skgeom._skgeom import properties


M = TypeVar('M', Mesh2, Mesh3)


class Mesh(Generic[M]):
    def __init__(self, mesh: M):
        self._mesh = mesh
        self.vertex_data = MeshData(mesh, properties.add_vertex_property, 'vertices')

    @staticmethod
    def from_polygon_soup(verts: A, faces: A, orient=True) -> Mesh:
        if verts.shape[1] == 2:
            mesh = Mesh2(verts, faces, orient)
        else:
            mesh = Mesh3(verts, faces, orient)
        return Mesh(mesh)

    @property
    def vertices(self):
        return array(self._mesh.vertices)

    @property
    def edges(self):
        return array(self._mesh.edges)

    @property
    def n_vertices(self) -> int:
        return self._mesh.n_vertices


class MeshData(Generic[M]):
    def __init__(self, mesh: M, add_fn, key_name: str):
        self._data = {}
        self._mesh = mesh
        self._add_fn = add_fn
        self._key_name = key_name

    def add_property(self, key: str, default: Any):
        pmap = self._add_fn(self._mesh, key, default)
        self._data[key] = pmap

    def __getitem__(self, item: str):
        return self._data[item]

    def __setitem__(self, key: str, value: ndarray):
        # Find or create the property map
        if key in self._data:
            pmap = self._data[key]
        else:
            default = zeros_like(value, shape=()).item()
            pmap = self._add_fn(self._mesh, key, default)
            self._data[key] = pmap

        # Set values
        keys = getattr(self._mesh, self._key_name)
        pmap[keys] = value
