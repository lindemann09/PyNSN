from re import A
from typing import List, Sequence
import numpy as np
from numpy.typing import NDArray
import pynsn
from pynsn._lib import np_coordinates, np_dots, np_rectangles, np_tools
from tests.shapes_test_picture import shapes_test_picture, Line

b = pynsn.Rectangle((55, 60), size=(40, 40), attribute="#FF0000")
a = pynsn.Dot((8, 9), diameter=88, attribute="#200800")

xy = np.array([[100, -1000], [100, 100],  [100, 100]])
sizes = np.array([[200, 200], [100, 100],  [100, 100]])
#xy = np.array([[50, 60]])
#sizes = np.array([[200, 200]])
dot_xy = np_tools.as_array2d_nrow(a.xy, n_rows=xy.shape[0])
dot_diameter = np_tools.as_array2d_nrow(900, n_rows=xy.shape[0])

over = np_rectangles.corner_inside_dot(rect_xy=xy, rect_sizes=sizes,
                                       dot_diameter=dot_diameter, dot_xy=dot_xy)
print(over)
exit()
# a_xy, a_size, b_xy = np_tools.all_as_array2d_equal_rows((a.xy, a.size, b.xy))

overlap = np_rectangles.corner_inside_dot(
    rect_xy=b.xy, rect_sizes=b.size,
    dot_xy=a.xy, dot_diameter=np.asarray([a.diameter]))
print(overlap)

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
shapes_test_picture((a, b, l))
