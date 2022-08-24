
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
aa = pynsn.Dot((55, 50), diameter=10, attribute="#200800")


b = pynsn.Rectangle((50, 50), size=(40, 40), attribute="#FF0000")
c = pynsn.Rectangle((-32, 60), size=(40, 40), attribute="#ccFF00")
d = pynsn.Rectangle((5, -71), size=(40, 40), attribute="#000088")


xy = np.array([b.xy, c.xy, d.xy])
sizes = np.array([b.size, c.size, d.size])
# xy = np.array([b.xy])
# sizes = np.array([b.size])
dot_xy = np_tools.as_array2d_nrow(a.xy, n_rows=xy.shape[0])
dot_diameter = np_tools.as_array2d_nrow(a.diameter, n_rows=xy.shape[0])


rel = spatial_relations.RectangleDotSpatRel(rect_xy=xy,
                                            rect_sizes=sizes,
                                            dot_xy=np.atleast_1d(a.xy),
                                            dot_diameter=np.atleast_1d(a.diameter))
rel.spatial_relations()

displace = rel.required_displacements(minimum_distance=1)
print(displace)
# FIXME required displacement not correct
b.xy = b.xy + displace[0, :]
c.xy = c.xy + displace[1, :]
d.xy = d.xy + displace[0, :]

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
shapes_test_picture((a, b, c, d))
