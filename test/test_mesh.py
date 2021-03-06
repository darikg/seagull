import pytest
import os

from tempfile import TemporaryDirectory
import numpy as np
from numpy import array, cos, sin,  pi, ones, arange, uint64
from numpy.testing import assert_array_equal

from skgeom.mesh import Mesh


def tetrahedron(scale=1.0, rot_z=0.0):
    verts = scale * array([[1, 1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, 1]], dtype='float')

    if rot_z:
        rot = array([[cos(rot_z), -sin(-rot_z), 0], [sin(rot_z), cos(-rot_z), 0], [0, 0, 1]])
        verts = verts @ rot.T

    faces = array([[2, 1, 0], [2, 3, 1], [3, 2, 0], [1, 3, 0]], dtype='int')

    return Mesh.from_polygon_soup(verts, faces, orient=True)
    # return verts, faces


@pytest.fixture
def mesh():
    return tetrahedron()


# @pytest.mark.parametrize('mesh', [tetrahedron(), triangle()])
# def test_property(mesh):
#     mesh.vertex_data['foo'] = arange(mesh.n_vertices)
#     verts = mesh.vertices
#     assert mesh.vertex_data['foo'][verts[0]] == 0
#     assert mesh.vertex_data['foo'][verts[1]] == 1
#
#     mesh.vertex_data['foo'] = 1 + arange(mesh.n_vertices)
#     assert mesh.vertex_data['foo'][verts[0]] == 1
#     assert mesh.vertex_data['foo'][verts[1]] == 2


# def test_corefine():
#     from skgeom._skgeom.mesh import Mesh3
#     m1 = tetrahedron()
#     nv0 = m1.n_vertices
#     m2 = tetrahedron(scale=0.9, rot_z=pi/3)
#     ecm = m1.edge_data.add_property('ecm', False)
#     m1.corefine(m2, dict(edge_is_constrained_map='ecm'))
#     nv1 = m1.n_vertices
#     assert nv1 > nv0
#     assert ecm[m1.edges].any()

