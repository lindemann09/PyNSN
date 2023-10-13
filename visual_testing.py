import pynsn
from tests.shapes_test_picture import shapes_test_picture, Line
import shapely
import numpy as np

d_big = pynsn.Dot((-8, 9), diameter=120, attribute="#00FF00")
da = pynsn.Dot((-20, 70), diameter=60, attribute="#222800")
db = pynsn.Dot((-70, -50), diameter=40, attribute="#000088")

ra = pynsn.Rectangle((120, -25), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-54, -45), size=(40, 40), attribute="#ccFF00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
objects = [d_big, ra, rb]


print([shapely.Point(o.xy) for o in objects])
polys = shapely.Polygon([o.xy for o in objects])
ch = polys.convex_hull

first = None
last = None
for xy in ch.exterior.coords:
    if first is None:
        first = xy
        last = xy
    else:
        l = Line(xy_a=last, xy_b=xy, colour="red")
        objects.append(l)
        last = xy
l = Line(xy_a=last, xy_b=first, colour="red")
objects.append(l)

shapes_test_picture(objects)

print(ch.centroid.x)
