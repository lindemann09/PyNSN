
import numpy as np
import pynsn
from pynsn._lib import np_tools, geometry
from pynsn import spatial_relations as sprel
from tests.shapes_test_picture import shapes_test_picture, Line

arr_xy = np.array([[10, 20], [100, 100], [99, 32]])
arr_sizes = np.array([[10, 2], [10, 100], [23, 4]])

arr_dia = np.array([10, 20])

da = pynsn.Dot((8, 9), diameter=120, attribute="#FF0000")
db = pynsn.Dot((70, 7), diameter=60, attribute="#222800")
dc = pynsn.Dot((-76, 5), diameter=40, attribute="#000088")


ra = pynsn.Rectangle((50, 50), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-32, 60), size=(40, 40), attribute="#ccFF00")
rc = pynsn.Rectangle((10, -40), size=(150, 50), attribute="#000088")

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")

shapes_test_picture((da, db, dc, l))


array_a = pynsn.DotArray()
array_a.add((db, dc))

array_b = pynsn.DotArray()
array_b.add((da))

sr = sprel.DotDot(array_a, array_b, reversed_relations=False)
print(f"dist: {sr.distances}\nangle {sr.angles}\nis_inside: {sr.is_inside()}\n")

x = np.array([sr.distances, sr.angles]).T
#print(f"{np.round(geometry.polar2cartesian(x), decimals=2)}")
