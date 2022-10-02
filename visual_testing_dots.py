import numpy as np
import pynsn
from pynsn._lib import np_tools, geometry
from pynsn import spatial_relations as sprel
from tests.shapes_test_picture import shapes_test_picture, Line

arr_xy = np.array([[10, 20], [100, 100], [99, 32]])
arr_sizes = np.array([[10, 2], [10, 100], [23, 4]])

arr_dia = np.array([10, 20])

d_big = pynsn.Dot((-8, 9), diameter=120, attribute="#FF0000")
da = pynsn.Dot((70, 70), diameter=60, attribute="#222800")
db = pynsn.Dot((-27, -50), diameter=40, attribute="#000088")


ra = pynsn.Rectangle((50, 50), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-50, -45), size=(40, 40), attribute="#ccFF00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
shapes_test_picture((d_big, da, db,  l))


# array_a = pynsn.RectangleArray()
# array_a.add((ra, rb))
# array_b = pynsn.RectangleArray()
# array_b.add(r_big)

array_a = pynsn.DotArray()
array_a.add((da, db))
array_b = pynsn.DotArray()
array_b.add(d_big)


sr = sprel.DotDot(array_a, array_b, A_relative_to_B=False)

print(f"dist: {sr.distances_rho}\nangle {np.rad2deg(sr.rho)}\nis_inside: {sr.is_inside()}\n"
      f"displ_dist: {sr.displacements_polar_old()}")

print("  ")
#print(f"{np.round(geometry.polar2cartesian(x), decimals=2)}")

displ = sr.displacements_cartesian(minimum_gap=0)
print(displ)
array_b.xy = array_b.xy + displ[1, :]

tmp = [l] + list(array_a.iter()) + list(array_b.iter())
shapes_test_picture(tmp, filename="shape_test2.png")

exit()
sr = sprel.DotDot(array_a, array_b, A_relative_to_B=True)
print(f"dist: {sr.distances_rho}\nangle {np.rad2deg(sr.rho)}\nis_inside: {sr.is_inside()}\n"
      f"displ_dist: {sr.displacements_polar_old()}")
