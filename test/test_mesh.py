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


def triangle(scale=1.0, rotate=0.0):
    verts = scale * array([[-1, -1], [1, -1], [0, 1]], dtype='float')

    if rotate:
        rot = array([[cos(rotate), -sin(-rotate)], [sin(rotate), cos(-rotate)]])
        verts = verts @ rot.T

    faces = array([[0, 1, 2]], dtype='int')

    return Mesh.from_polygon_soup(verts, faces, orient=True)


@pytest.fixture
def mesh():
    return tetrahedron()


@pytest.fixture
def mesh2():
    return triangle()


@pytest.mark.parametrize('mesh', [tetrahedron(), triangle()])
def test_property(mesh):
    mesh.vertex_data['foo'] = arange(mesh.n_vertices)
    verts = mesh.vertices
    assert mesh.vertex_data['foo'][verts[0]] == 0
    assert mesh.vertex_data['foo'][verts[1]] == 1

    mesh.vertex_data['foo'] = 1 + arange(mesh.n_vertices)
    assert mesh.vertex_data['foo'][verts[0]] == 1
    assert mesh.vertex_data['foo'][verts[1]] == 2
