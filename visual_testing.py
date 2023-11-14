from pynsn import random
from pynsn._misc import cartesian2polar
from pynsn.random._rng import WalkAround
import time
from timeit import default_timer as timer
import pynsn
from pynsn.image import pil_image, svg_file
import numpy as np

from pynsn import random




da = pynsn.Dot((-20, 120), diameter=5,  attribute=pynsn.Colour("black"))
db = pynsn.Ellipse((20, 17), size=(120, 50), attribute="#FF0000")
dc = pynsn.Ellipse((120, 57), size=(60, 120), attribute="#00F000")  # big dot

ra = pynsn.Rectangle(np.array((-70, -175)), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-50, -45), size=(10, 10), attribute="#cc0F00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

rnd_ell = random.RndEllipse(width=(40.8, 10.4), height=(20, 50),
                            attributes=["green", "black", "orange", "red"])
rnd_ell = random.RndDot(diameter=(40.8, 10),
                            attributes=["green", "black", "orange", "red"])

nsn = pynsn.NSNStimulus(
    #target_area_shape=pynsn.Dot((0, 0), diameter=500, attribute="#00FFFF"),
    target_area_shape=pynsn.Rectangle((0, 0), size=(400, 500)),
    min_distance_target_area=10,
    min_distance=2)


# random dot
nsn.add_shape(rnd_ell, n=20, random_position=True)

#nsn.add_shape_list([db,  dc,  ra], ignore_overlaps=False)
#nsn.fix_overlaps()

print(nsn.properties_txt(short_format=True))
print(nsn.contains_overlaps())
print(nsn.properties.numerosity)

nsn.sort_by_excentricity()
a = pil_image.create(nsn)
a.save("shapes_test.png")

d = nsn.tojson(filename="demo.json", tabular=True)