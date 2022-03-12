import pytest
import os

from tempfile import TemporaryDirectory
import numpy as np
from numpy import array, cos, sin,  pi, ones, arange, uint64
from numpy.testing import assert_array_equal

from skgeom._skgeom.surface_mesh import Mesh


def tetrahedron(scale=1.0, rot_z=0.0):
    verts = scale * array([[1, 1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, 1]], dtype='float')

    if rot_z:
        rot = array([[cos(rot_z), -sin(-rot_z), 0], [sin(rot_z), cos(-rot_z), 0], [0, 0, 1]])
        verts = verts @ rot

    faces = array([[2, 1, 0], [2, 3, 1], [3, 2, 0], [1, 3, 0]], dtype='int')

    return verts, faces


@pytest.fixture
def mesh():
    return Mesh(*tetrahedron(), True)


def test_scalar_property(mesh):
    v = mesh.vertices[0]
    pmap = mesh.add_vertex_property('foo', -1)
    assert pmap[v] == -1
    pmap[v] = 1
    assert pmap[v] == 1
    assert pmap[mesh.vertices[1]] == -1


def test_vector_property(mesh):
    verts = mesh.vertices
    pmap = mesh.add_vertex_property('foo', 1)
    assert_array_equal(pmap[verts], ones(mesh.n_vertices))

    vals = arange(mesh.n_vertices)
    pmap[verts] = vals
    assert_array_equal(pmap[verts], vals)


@pytest.mark.parametrize('ext', ('off', 'ply'))
def test_io(mesh, ext):
    with TemporaryDirectory() as d:
        file = os.path.join(d, 'mesh.' + ext)
        getattr(mesh, 'write_' + ext)(file)
        _ = Mesh(file)


def test_corefine():
    m0 = Mesh(*tetrahedron(), True)
    m1 = Mesh(*tetrahedron(scale=0.9,  rot_z=pi/3), True)
    m0.corefine(m1)


def test_corefine_tracked():
    m1 = Mesh(*tetrahedron(), True)
    vert_ids1 = m1.add_vertex_property('corefined_idx', -1)
    ecm1 = m1.add_edge_property('constrained', False)

    m2 = Mesh(*tetrahedron(scale=0.9,  rot_z=pi/3), True)
    vert_ids2 = m2.add_vertex_property('corefined_idx', -1)
    ecm2 = m2.add_edge_property('constrained', False)
    m1.corefine(vert_ids1, ecm1, m2, vert_ids2, ecm2)

    # Corefined vert indices should match
    v1, v2 = m1.vertices, m2.vertices
    ci1, ci2 = vert_ids1[v1], vert_ids2[v2]
    assert (ci1 > -1).any()
    assert (ci2 > -1).any()
    assert_array_equal(sorted(ci1[ci1 > -1]), sorted(ci2[ci2 > -1]))

    # The intersection edges should've been marked as constrained
    for (m, ecm) in zip((m1, m2), (ecm1, ecm2)):
        assert ecm[m.edges].any()


@pytest.mark.parametrize('op', ['difference', 'union', 'intersection'])
def test_binary_op(op):
    m0 = Mesh(*tetrahedron(), True)
    m1 = Mesh(*tetrahedron(scale=0.9, rot_z=pi / 3), True)
    result = getattr(m0, op)(m1)
    assert isinstance(result, Mesh)


def test_remesh_simple(mesh):
    nv0 = mesh.n_vertices
    mesh.remesh(mesh.faces[:2], 0.5, 1)
    assert mesh.n_vertices > nv0


def test_remesh_tracked(mesh):
    nv0 = mesh.n_vertices
    touched = mesh.add_vertex_property('touched', False)
    mesh.remesh(mesh.faces[:2], 0.5, 1, touched)
    assert mesh.n_vertices > nv0
    is_touched = touched[mesh.vertices]
    assert is_touched.any()
    assert ~is_touched.all()


def test_fair(mesh):
    mesh.fair(mesh.vertices[:2], 0)

