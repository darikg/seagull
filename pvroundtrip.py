from pyvista import Sphere
from seagullmesh.mesh import Mesh3

s = Sphere().clean().triangulate()
s2 = Mesh3.from_pyvista(s).to_pyvista()
s2.plot()