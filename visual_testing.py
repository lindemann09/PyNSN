
from typing import List, Sequence
import numpy as np
from numpy.typing import NDArray
import pynsn
from pynsn._lib import geometry, np_tools, spatial_relations, matrix_spatial_relations
from tests.shapes_test_picture import shapes_test_picture, Line

arr_xy = np.array([[10, 20], [100, 100], [99, 32]])
arr_sizes = np.array([[10, 2], [10, 100], [23, 4]])

points = [[1, 10], [37, 1], [2.5, 2], [2.5, 1]]


def p4(p1, p2, p3):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    dx, dy = x2-x1, y2-y1
    det = dx*dx + dy*dy
    a = (dy*(y3-y1)+dx*(x3-x1))/det

    return x1+a*dx, y1+a*dy


# print(np.ones((3, 1)) * np.atleast_2d(points[0]))
#print(p4(points[0], points[1], points[2]))
p4 = line_point_othogonal(points[0], points[1], points[2])
dxy = np.asarray(points[2]) - p4
print(np.hypot(dxy[0], dxy[1]))

exit()


arr_dia = np.array([10, 20])
a = pynsn.Dot((8, 9), diameter=120, attribute="#200800")
aa = pynsn.Dot((55, 50), diameter=10, attribute="#200800")


b = pynsn.Rectangle((50, 50), size=(40, 40), attribute="#FF0000")
c = pynsn.Rectangle((-32, 60), size=(40, 40), attribute="#ccFF00")
d = pynsn.Rectangle((5, -60), size=(150, 40), attribute="#000088")


xy = np.array([b.xy, c.xy, d.xy])
sizes = np.array([b.size, c.size, d.size])
# xy = np.array([b.xy])
# sizes = np.array([b.size])
dot_xy = np_tools.as_array2d_nrow(a.xy, n_rows=xy.shape[0])
dot_diameter = np_tools.as_array2d_nrow(a.diameter, n_rows=xy.shape[0])


rel = spatial_relations.RectangleDotSpatRel(rect_xy=xy,
                                            rect_sizes=sizes,
                                            dot_xy=np.atleast_2d(
                                                (a.xy, a.xy, a.xy)),
                                            dot_diameter=np.atleast_1d(
                                                (a.diameter, a.diameter, a.diameter)))
print(rel.spatial_relations())
print(rel._cardinal_point_edge_relations())


displace = rel.required_displacements(minimum_distance=1)
print(displace)
# FIXME required displacement not correct
# b.xy = b.xy + displace[0, :]
# c.xy = c.xy + displace[1, :]
# d.xy = d.xy + displace[0, :]

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
shapes_test_picture((a, b, c, d))
