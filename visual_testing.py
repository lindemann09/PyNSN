import time
from timeit import default_timer as timer
import pynsn
from pynsn.image import pil_image
import shapely
from tests.shapes_test_picture import shapes_test_picture, Line
import numpy as np

d_big = pynsn.Dot((-8, -67), diameter=120, attribute="#00FF00")
da = pynsn.Dot((-20, 120), diameter=60, attribute="#222800")
db = pynsn.Dot((-70, -50), diameter=10, attribute="#000088")

ra = pynsn.Rectangle(np.array((120, -175)), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-54, -45), size=(10, 10), attribute="#cc0F00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

nsn = pynsn.NSNStimulus(
    target_area=pynsn.Dot((0, 0), diameter=500, attribute="#00FFFF")
)

start = timer()
nsn.add_somewhere(db, n=400, ignore_overlaps=False)
# nsn.add([da, db, ra, rb])
# print(nsn.properties_txt(extended_format=True))
end = timer()
print(f"adding :{end - start}")

col = pil_image.ImageColours(field_area="red")
a = pil_image.create(nsn, colours=col)
a.save("shapes_test2.png")

# print(nsn.target_area.polygon.covers(nsn.polygons))
# print(pynsn.spatial_relations.dwithin(nsn, dist=nsn.min_distance))
