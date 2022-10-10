
import numpy as np
import pynsn
from pynsn._lib import np_tools, geometry
from pynsn import spatial_relations as sprel
from tests.shapes_test_picture import shapes_test_picture, Line

arr_xy = np.array([[10, 20], [100, 100], [99, 32]])
arr_sizes = np.array([[10, 2], [10, 100], [23, 4]])

arr_dia = np.array([10, 20])

d_big = pynsn.Dot((-8, 90), diameter=120, attribute="#FF0000")
da = pynsn.Dot((35, -40), diameter=60, attribute="#222800")
db = pynsn.Dot((-120, 55), diameter=40, attribute="#00FFF0")


ra = pynsn.Rectangle((50, 50), size=(40, 40), attribute="#FF000")
rb = pynsn.Rectangle((-10, 45), size=(40, 40), attribute="#00FF00")
r_big = pynsn.Rectangle((0, 0), size=(150, 60), attribute="#000088")

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")

a_relative_to_b = False
array_rect = pynsn.RectangleArray()
array_dot = pynsn.DotArray()

if a_relative_to_b:
    # A_rel_b
    array_rect.add((r_big, rb))
    array_dot.add(d_big)
else:
    array_rect.add(r_big)
    array_dot.add((da, db))

a = list(array_rect.iter()) + list(array_dot.iter())
shapes_test_picture(a)

sr = sprel.RectangleDot(array_rect, array_dot, a_relative_to_b=a_relative_to_b)
print(f"dist: {sr.distances_rho}\nangle {sr.rho}\nis_inside: {sr.is_inside()}\n"
      f"displ_dist: {sr._spread_distances_rho()}")

print("  ")
# print(f"{np.round(geometry.polar2cartesian(x), decimals=2)}")


displ = sr.spread(displ_type=sprel.SpreadTypes.SHORTEST,
                  minimum_gap=0)
print(displ)
array_dot.xy = array_dot.xy + displ

tmp = []
tmp.extend(array_rect.iter())
tmp.extend(array_dot.iter())
shapes_test_picture(tmp, filename="shape_test2.png")


#sr = sprel.RectangleDot(array_rect, array_dot, a_relative_to_b=True)
# print(f"dist: {sr.distances_rho}\nangle {sr.rho}\nis_inside: {sr.is_inside()}\n"
#      f"displ_dist: {sr.displacement_distances_rho()}")
