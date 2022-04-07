import pytest
import os

from tempfile import TemporaryDirectory, TemporaryFile
from pathlib import Path

import numpy as np
from numpy import array, cos, sin,  pi, ones, arange, uint64, full, zeros
from numpy.testing import assert_array_equal

from seagullmesh import Mesh3, MeshData


def tetrahedron(scale=1.0, rot_z=0.0):
    verts = scale * array([[1, 1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, 1]], dtype='float')

    if rot_z:
        rot = array([[cos(rot_z), -sin(-rot_z), 0], [sin(rot_z), cos(-rot_z), 0], [0, 0, 1]])
        verts = verts @ rot.T

    faces = array([[2, 1, 0], [2, 3, 1], [3, 2, 0], [1, 3, 0]], dtype='int')
    return verts, faces


def test_from_polygon_soup():
    verts, faces = tetrahedron()
    mesh = Mesh3.from_polygon_soup(verts, faces)
    assert mesh.n_vertices == 4 and mesh.n_faces == 4


@pytest.mark.parametrize('file', ['armadillo.off', 'sphere.ply'])
def test_from_ply(file):
    file = Path(__file__).parent / 'assets' / file
    assert file.exists()
    mesh = Mesh3.from_file(str(file))


@pytest.mark.parametrize('ext', ['ply', 'off'])
def test_to_file(ext):
    mesh = Mesh3.from_polygon_soup(*tetrahedron())
    with TemporaryDirectory() as d:
        file = str(Path(d) / f'mesh.{ext}')
        mesh.to_file(file)


def test_pyvista_roundtrip():
    from pyvista import Sphere
    pvmesh0 = Sphere().clean().triangulate()
    mesh = Mesh3.from_pyvista(pvmesh0)
    pvmesh1 = mesh.to_pyvista()


@pytest.mark.parametrize('key_type', ['vertex', 'face', 'edge', 'halfedge'])
@pytest.mark.parametrize('val_type', [int, bool])
def test_scalar_properties(key_type, val_type):
    mesh = Mesh3.from_polygon_soup(*tetrahedron())
    d: MeshData = getattr(mesh, key_type + '_data')

    d['foo'] = full(d.n_mesh_keys, val_type(0))

    keys = d.mesh_keys
    key = keys[0]
    d['foo'][key] = val_type(1)
    val = d['foo'][key]
    assert val == val_type(1)

    d['foo'][keys[:2]] = [val_type(1), val_type(1)]
    assert d['foo'][keys[0]] == val_type(1) and d['foo'][keys[1]] == val_type(1)


def corefine_meshes():
    m1 = Mesh3.from_polygon_soup(*tetrahedron())
    m2 = Mesh3.from_polygon_soup(*tetrahedron(scale=0.9, rot_z=pi/3))
    return m1, m2


def test_corefine():
    m1, m2 = corefine_meshes()
    nv1_orig, nv2_orig = m1.n_vertices, m2.n_vertices
    m1.corefine(m2)

    nv1, nv2 = m1.n_vertices, m2.n_vertices
    assert nv1 > nv1_orig and nv2 > nv2_orig


@pytest.mark.parametrize('op', ['union', 'intersection', 'difference'])
@pytest.mark.parametrize('inplace', [False, True])
def test_boolean_ops(op, inplace):
    m1, m2 = corefine_meshes()
    m3 = getattr(m1, op)(m2, inplace=inplace)
