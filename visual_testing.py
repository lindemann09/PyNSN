
from typing import List, Sequence
import numpy as np
from numpy.typing import NDArray
import pynsn
from pynsn._lib import geometry, np_tools, spatial_relations, matrix_spatial_relations
from tests.shapes_test_picture import shapes_test_picture, Line

arr_xy = np.array([[10, 20], [100, 100], [99, 32]])
arr_sizes = np.array([[10, 2], [10, 100], [23, 4]])


arr_dia = np.array([10, 20])
a = pynsn.Dot((8, 9), diameter=120, attribute="#200800")

b = pynsn.Rectangle((50, 50), size=(40, 40), attribute="#FF0000")
c = pynsn.Rectangle((32, 60), size=(40, 40), attribute="#ccFF00")
d = pynsn.Rectangle((5, -71), size=(40, 40), attribute="#000088")

xy = np.array([b.xy, c.xy, d.xy])
sizes = np.array([b.size, c.size, d.size])
# xy = np.array([b.xy])
# sizes = np.array([b.size])
dot_xy = np_tools.as_array2d_nrow(a.xy, n_rows=xy.shape[0])
dot_diameter = np_tools.as_array2d_nrow(a.diameter, n_rows=xy.shape[0])


mtx = matrix_spatial_relations.RectangleMatrixSpatRel(xy, sizes)
displace = mtx.required_displacements(minimum_distance=2)
print(displace)

# FIXME required displacement not correct
c.xy = c.xy + np.ceil(displace[0, (2, 3)])


print(c.distance(b))

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
shapes_test_picture((a, b, c, d))
