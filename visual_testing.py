from re import A
from typing import List, Sequence
import numpy as np
from numpy.typing import NDArray
import pynsn
from pynsn._lib import np_coordinates, np_dots, np_rectangles, np_tools
from tests.shapes_test_picture import shapes_test_picture, Line

arr_xy = np.array([[10, 20], [100, 100], [99, 32]])
arr_sizes = np.array([[10, 2], [10, 100], [23, 4]])


arr_dia = np.array([10, 20])
d = pynsn.Rectangle((5, -71), size=(40, 40), attribute="#FF0088")
c = pynsn.Rectangle((-65, 60), size=(40, 40), attribute="#cc0011")
b = pynsn.Rectangle((50, 50), size=(40, 40), attribute="#FF0000")
a = pynsn.Dot((8, 9), diameter=120, attribute="#200800")

xy = np.array([b.xy, c.xy, d.xy])
sizes = np.array([b.size, c.size, d.size])
#xy = np.array([b.xy])
#sizes = np.array([b.size])
dot_xy = np_tools.as_array2d_nrow(a.xy, n_rows=xy.shape[0])
dot_diameter = np_tools.as_array2d_nrow(a.diameter, n_rows=xy.shape[0])

rdr = np_rectangles.RectangleDotRelations(xy, sizes, dot_xy, dot_diameter)

print(rdr.overlap(minimum_distance=2))


l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
shapes_test_picture((a, b, c, d))
