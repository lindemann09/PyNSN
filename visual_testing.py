
import numpy as np
import pynsn
from pynsn._lib import np_tools, geometry
from pynsn import spatial_relations as sprel
from tests.shapes_test_picture import shapes_test_picture, Line

arr_xy = np.array([[10, 20], [100, 100], [99, 32]])
arr_sizes = np.array([[10, 2], [10, 100], [23, 4]])

arr_dia = np.array([10, 20])

da = pynsn.Dot((-8, 9), diameter=120, attribute="#FF0000")
db = pynsn.Dot((70, 7), diameter=60, attribute="#222800")
dc = pynsn.Dot((-76, 5), diameter=40, attribute="#000088")


ra = pynsn.Rectangle((50, 50), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-65, -45), size=(40, 40), attribute="#ccFF00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
shapes_test_picture((r_big, ra, rb,  l))


array_a = pynsn.RectangleArray()
array_a.add((ra, rb))

array_b = pynsn.RectangleArray()
array_b.add(r_big)

sr = sprel.RectangleRectangle(array_a, array_b, reversed_relations=True)
print(f"dist: {sr.distances}\nangle {sr.angle}\nis_inside: {sr.is_inside()}\n"
      f"displ_dist: {sr.displacement_distances()}")

print(np.rad2deg(sr.angle))
#print(f"{np.round(geometry.polar2cartesian(x), decimals=2)}")

displ = sr.required_displacements(minimum_distance=0)
array_a.xy = array_a.xy + displ

tmp = [r_big]
tmp.extend(array_a.iter())
shapes_test_picture(tmp, filename="shape_test2.png")


sr = sprel.RectangleRectangle(array_a, array_b, reversed_relations=True)
print(f"dist: {sr.distances}\nangle {sr.angle}\nis_inside: {sr.is_inside()}\n"
      f"displ_dist: {sr.displacement_distances()}")
