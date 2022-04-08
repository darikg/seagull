import pytest
import os

from tempfile import TemporaryDirectory
import numpy as np
from numpy import array, cos, sin,  pi, ones, arange, uint64, zeros
from numpy.testing import assert_array_equal

from skgeom._skgeom.surface_mesh import Mesh3


def tetrahedron(scale=1.0, rot_z=0.0):
    verts = scale * array([[1, 1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, 1]], dtype='float')

    if rot_z:
        rot = array([[cos(rot_z), -sin(-rot_z), 0], [sin(rot_z), cos(-rot_z), 0], [0, 0, 1]])
        verts = verts @ rot

    faces = array([[2, 1, 0], [2, 3, 1], [3, 2, 0], [1, 3, 0]], dtype='int')

    return verts, faces


# m1 = Mesh(*tetrahedron(), True)
# vert_ids1 = m1.add_vertex_property('corefined_idx', -1)
# ecm1 = m1.add_edge_property('constrained', False)
#
# m2 = Mesh(*tetrahedron(scale=0.9, rot_z=pi / 3), True)
# vert_ids2 = m2.add_vertex_property('corefined_idx', -1)
# ecm2 = m2.add_edge_property('constrained', False)
# m1.corefine(vert_ids1, ecm1, m2, vert_ids2, ecm2)
#
# e1, e2 = array(m1.edges), array(m2.edges)
# ci1 = ecm1[e1]
# print(m1.edge_vertices(e1[ci1]))


def make_mesh():
    m = Mesh3(*tetrahedron(), True)
    m.add_edge_property('corefined_idx', False)
    return m


mesh = make_mesh()
print(mesh.expand_selection(mesh.vertices[:2]))
# tree = mesh.aabb_tree()
# face_idx, bary_coords = mesh.locate_with_aabb_tree(tree, zeros((1, 3)))
# # mesh = Mesh(*tetrahedron(), True)
# # pmap = mesh.add_edge_property('corefined_idx', False)
# pmap = mesh.get_edge_property('corefined_idx')
# print(pmap[mesh.edges])
