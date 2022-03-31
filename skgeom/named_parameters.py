
# def vertex_point_map(val)
from skgeom.mesh import Mesh, properties


def edge_is_constrained_map(mesh: Mesh, val):
    if isinstance(val, str):
        return mesh.edge_data.get_or_create_property(val, default=False)
    else:
        return val


