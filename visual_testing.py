from pynsn import random
import pynsn
from pynsn.image import pil_image
import numpy as np


da = pynsn.Dot((-20, 120), diameter=5,  attribute=pynsn.Colour("black"))
db = pynsn.Ellipse((20, 17), size=(120, 50), attribute="#FF0000")
dc = pynsn.Ellipse((120, 57), size=(60, 120), attribute="#00F000")  # big dot

ra = pynsn.Rectangle(np.array((-70, -175)), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-50, -45), size=(10, 10), attribute="#cc0F00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

rnd_ell = random.RndRectangle(width=(40.8, 10.4), height=(20, 50),
                              attributes=["green", "black", "orange", "red"])
rnd_dot = random.RndDot(diameter=(40.8, 10),
                        attributes=["green", "black", "orange", "red"])

if False:
    nsn = pynsn.NSNStimulus(
        # target_area_shape=pynsn.Dot((0, 0), diameter=500, attribute="#00FFFF"),
        target_area_shape=pynsn.Rectangle((0, 0), size=(400, 500)),
        min_distance_target_area=10,
        min_distance=2)
    # random dot
    nsn.add_shapes(rnd_ell, n=20, random_position=True)
    print(nsn.properties_txt(short_format=True))
    print(nsn.contains_overlaps())
    print(nsn.properties.numerosity)

if True:
    factory = pynsn.StimulusFactory(
        target_area_shape=pynsn.Rectangle((0, 0), size=(500, 400)),
        min_distance_target_area=10,
        min_distance=2)
    factory.add(rnd_ell, 10)
    factory.add(rnd_dot, 10)

    nsn = factory.get()
    factory.tojson(filename="demo.json")

a = pil_image.create(nsn)
a.save("shapes_test2.png")


sizes = nsn.get_sizes()
sizes[5:9, :] = sizes[5:9, :] * 3
nsn.set_sizes(sizes)
print(nsn.get_sizes())

a = pil_image.create(nsn)
a.save("shapes_test.png")
