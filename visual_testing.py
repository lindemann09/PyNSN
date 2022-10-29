
import numpy as np
import pynsn
from pynsn._lib import np_tools, geometry
from pynsn import spatial_relations as sprel
from tests.shapes_test_picture import shapes_test_picture, Line

arr_xy = np.array([[10, 20], [100, 100], [99, 32]])
arr_sizes = np.array([[10, 2], [10, 100], [23, 4]])

arr_dia = np.array([10, 20])

d_big = pynsn.Dot((-8, 90), diameter=120, attribute="#FF00FF")
da = pynsn.Dot((35, -40), diameter=60, attribute="#002800")
db = pynsn.Dot((-120, 55), diameter=40, attribute="#00FFF0")

ra = pynsn.Rectangle((50, 50), size=(40, 40), attribute="#0000FF")
rb = pynsn.Rectangle((-10, 45), size=(40, 40), attribute="#00FF00")
rc = pynsn.Rectangle((-75, -32), size=(40, 40), attribute="#0000FF")
rd = pynsn.Rectangle((10, -45), size=(40, 40), attribute="#00FF00")
r_big = pynsn.Rectangle((0, 0), size=(150, 60), attribute="#000FF0")

a_relative_to_b = True
array_rect = pynsn.RectangleArray()
array_dot = pynsn.DotArray()

if a_relative_to_b:
    # A_rel_b
    array_rect.add((ra, rb, rc, rd))
    array_dot.add(d_big)
else:
    array_rect.add(r_big)
    array_dot.add((da, db))

a = list(array_rect.iter()) + list(array_dot.iter())
shapes_test_picture(a, reverse_order=True)

sr = sprel.relations(array_rect, array_dot, a_relative_to_b=a_relative_to_b)
print(type(sr))

print(f"dist: {sr.distances_rho}\nangle {sr.rho}\nis_inside: {sr.is_inside()}\n")
print("  ")

# print(f"{np.round(geometry.polar2cartesian(x), decimals=2)}")

displ = sr.gather(minimum_gap=10)
print("---")
print(displ)
array_rect.xy = array_rect.xy + displ

tmp = []
tmp.extend(array_rect.iter())
tmp.extend(array_dot.iter())
shapes_test_picture(tmp, filename="shape_test2.png", reverse_order=True)


#sr = sprel.RectangleDot(array_rect, array_dot, a_relative_to_b=True)
# print(f"dist: {sr.distances_rho}\nangle {sr.rho}\nis_inside: {sr.is_inside()}\n"
#      f"displ_dist: {sr.displacement_distances_rho()}")
