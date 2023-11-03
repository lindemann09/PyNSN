import time
from timeit import default_timer as timer
import pynsn
from pynsn.image import pil_image
import numpy as np



da = pynsn.Dot((-20, 120), diameter=10, attribute="#222800")
db = pynsn.Dot((20, 17), diameter=60, attribute="#FF0000")
dc = pynsn.Ellipse((120, 57), size=(60, 120), attribute="#00F000")  # big dot

ra = pynsn.Rectangle(np.array((-70, -175)), size=(40, 40), attribute="#FF0000")
rb = pynsn.Rectangle((-50, -45), size=(10, 10), attribute="#cc0F00")
r_big = pynsn.Rectangle((10, -40), size=(150, 60), attribute="#000088")

nsn = pynsn.NSNStimulus(
    # target_area=pynsn.Dot((0, 0), diameter=500, attribute="#00FFFF")
    target_area=pynsn.Ellipse((0, 0), size=(400, 500), attribute="#00FFFF"),
    min_distance_target_area=30,
    min_distance=2
)
nsn.add([db,  dc, ra])


if True:
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
    start = timer()
    nsn.add_somewhere(da, n=200, ignore_overlaps=False)
    end = timer()
    print(f"adding :{end - start}")
    print(nsn.properties_txt(extended_format=True))


col = pil_image.ImageColours(field_area="red")
a = pil_image.create(nsn, colours=col)
a.save("shapes_test2.png")

