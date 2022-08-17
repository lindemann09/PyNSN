from typing import List, Sequence
import numpy as np
from numpy.typing import NDArray
import pynsn
from pynsn._lib import np_coordinates, np_dots, np_rectangles, np_tools
from tests.shapes_test_picture import shapes_test_picture, Line

b = pynsn.Rectangle((55, 60), size=(40, 40), attribute="#FF0000")
a = pynsn.Dot((00, 00), diameter=88, attribute="#200800")

# a_xy, a_size, b_xy = np_tools.all_as_array2d_equal_rows((a.xy, a.size, b.xy))

overlap = np_rectangles.overlap_with_dots(
    rect_xy=b.xy, rect_sizes=b.size,
    dot_xy=a.xy, dot_diameter=np.asarray([a.diameter]))
print(overlap)

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
shapes_test_picture((a, b, l))
