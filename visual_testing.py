import pynsn
from pynsn.image import pil_image

from tests.shapes_test_picture import shapes_test_picture, Line
import numpy as np

d_big = pynsn.Dot((-8, -67), diameter=120, attribute="#00FF00")
da = pynsn.Dot((-20, 70), diameter=60, attribute="#222800")
db = pynsn.Dot((-70, -50), diameter=40, attribute="#000088")

ra = pynsn.Rectangle(np.array((120, -25)), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-54, -45), size=(40, 40), attribute="#ccFF00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

nsn = pynsn.NSNStimulus(
    target_area=pynsn.Dot((0, 0), diameter=500, attribute="#00FFFF")
)

nsn.shapes.add([da, db, ra, rb])
col = pil_image.ImageColours(field_area="red")
a = pil_image.create(nsn, colours=col)
a.save("shapes_test2.png")


l = Line(xy_a=(200, 200), xy_b=(-300, 0), colour="red")
objects = [d_big, ra, rb]

shapes_test_picture(objects)
