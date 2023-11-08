from pynsn import random
from pynsn._misc import cartesian2polar
from pynsn.random._rng import WalkAround
import time
from timeit import default_timer as timer
import pynsn
from pynsn.image import pil_image, svg_file
import numpy as np


da = pynsn.Dot((-20, 120), diameter=10)
db = pynsn.Ellipse((20, 17), size=(120, 50), attribute="#FF0000")
dc = pynsn.Ellipse((120, 57), size=(60, 120), attribute="#00F000")  # big dot

ra = pynsn.Rectangle(np.array((-70, -175)), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-50, -45), size=(10, 10), attribute="#cc0F00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

nsn = pynsn.NSNStimulus(
    # target_area=pynsn.Dot((0, 0), diameter=500, attribute="#00FFFF")
    target_area_shape=pynsn.Rectangle(
        (0, 0), size=(400, 500), attribute="#00FFFF"),
    min_distance_target_area=10,
    min_distance=2
)

if False:
    dist = 2
    a = db
    b = dc
    print(nsn.matrix_distances())

    # sp_result = shapely.distance(a.polygon, b.polygon)
    # geo_result = sgeo.distance(a, b)
    # arr_result = sgeo.distance_circ_dot_array(a,
    #                                           dots_xy=np.atleast_2d(b.xy),
    #                                           dots_diameter=np.atleast_2d(b.diameter))
    # arr_result2 = sgeo.distance_circ_ellipse_array(a,
    #                                                ellipses_xy=np.atleast_2d(
    #                                                    b.xy),
    #                                                ellipse_sizes=np.atleast_2d(b.size))
    # print((sp_result, geo_result, arr_result, arr_result2))

# random dot
if True:
    nsn.add_somewhere(da, n=20, ignore_overlaps=False)
nsn.add([db,  dc, ra])
# nsn.colours.object_default = "red"
# nsn.colours.convex_hull = "gray"

nsn.sort_by_excentricity()
a = pil_image.create(nsn)
a.save("shapes_test.png")
