

import numpy as np
import pynsn
from pynsn._lib import np_tools, geometry
from pynsn import spatial_relations as sprel
from tests.shapes_test_picture import shapes_test_picture, Line

arr_xy = np.array([[10, 20], [100, 100], [99, 32]])
arr_sizes = np.array([[10, 2], [10, 100], [23, 4]])


d_big = pynsn.Dot((0, 0), diameter=220, attribute="black")
yd = pynsn.Dot((50, 50), diameter=4, attribute="red")
shapes_test_picture([d_big, yd], reverse_order=False)

xy = np.atleast_2d(yd.xy)
rho = np.deg2rad(230.00)
intersection = geometry.line_circle_intersection(
    line_points=xy,
    line_directions=np.asarray([rho]),
    circle_center=np.asarray([d_big.xy]),                                    circle_radii=np.asarray([d_big.diameter])/2)

print(intersection)


displ = intersection - xy
yd.xy = yd.xy + displ[0]

shapes_test_picture([d_big, yd], reverse_order=False,
                    filename="shape_test2.png")

exit()


arr_dia = np.array([10, 20])

d_big = pynsn.Dot((-10, 90), diameter=220, attribute="black")
da = pynsn.Dot((35, -40), diameter=60, attribute="#002800")
db = pynsn.Dot((-120, 55), diameter=40, attribute="#00FFF0")
dc = pynsn.Dot((20, -55), diameter=20, attribute="#FFFFF0")

ra = pynsn.Rectangle((50, 50), size=(50, 40), attribute="blue")
rb = pynsn.Rectangle((-30, 25), size=(90, 90), attribute="green")
rc = pynsn.Rectangle((-75, -32), size=(40, 40), attribute="yellow")
rd = pynsn.Rectangle((0, -45), size=(100, 40), attribute="magenta")
r_big = pynsn.Rectangle((0, 0), size=(150, 60), attribute="#000FF0")

a_relative_to_b = True
array_rect = pynsn.RectangleArray()
array_dot = pynsn.DotArray()

if a_relative_to_b:
    # A_rel_b
    array_rect.add((ra, rb, rc, rd))
    # array_rect.add((ra, rb))
    array_dot.add(d_big)
else:
    array_rect.add(r_big)
    array_dot.add((da, db, dc))

a = list(array_rect.iter()) + list(array_dot.iter())
shapes_test_picture(a, reverse_order=True)


sr = sprel.relations(a_array=array_rect, b_array=array_dot,
                     a_relative_to_b=True)
# print(f"{np.round(geometry.polar2cartesian(x), decimals=2)}")
displ = sr.spread(polar=False, minimum_gap=1, cardinal_axis_only=False)
# displ = sr.gather(polar=False, minimum_gap=1)
print(displ)
# displ = sr.gather(minimum_gap=0)
# print("---")
# print(displ)
array_rect.xy = array_rect.xy + displ

tmp = []
tmp.extend(array_rect.iter())
tmp.extend(array_dot.iter())
shapes_test_picture(tmp, filename="shape_test2.png", reverse_order=True)


# sr = sprel.RectangleDot(array_rect, array_dot, a_relative_to_b=True)
# print(f"dist: {sr.distances_rho}\nangle {sr.rho}\nis_inside: {sr.is_inside()}\n"
#      f"displ_dist: {sr.displacement_distances_rho()}")
