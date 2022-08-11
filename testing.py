"""_summary_"""

import time
from matplotlib import pyplot as pp
import numpy as np
from numpy.typing import NDArray, ArrayLike

import pynsn
import pynsn as nsn
from pynsn import distributions as distr
from pynsn._lib.misc import numpy_array_2d
from pynsn._shapes.rectangle import Rectangle
from pynsn.image import svg_file, pil_image, mpl_figure
from pynsn._lib import geometry

seed = 921
nsn.init_random_generator(seed)

my_colours = nsn.ImageColours(  # target_area="#EEEEEE",
    background="",
    opacity_object=0.9,
    default_object_colour="darkmagenta",
    # field_area_positions="magenta",
    # field_area="gray",
    # center_of_positions="red",
    # center_of_mass="magenta"
)


factory = nsn.NSNFactory(min_distance_between_objects=2,
                         #target_area=pynsn.Rectangle(size=(200, 400)),
                         target_area=pynsn.Dot(diameter=400),
                         min_distance_area_boarder=2)

factory.set_appearance_dots(diameter=(2, 3, 4),
                            attributes=distr.Levels(["blue", "green"],
                                                    exact_weighting=True))

factory.set_appearance_rectangles(width=(2, 3, 4), proportion=1,
                                  attributes=distr.Levels(["blue", "green"],
                                                          exact_weighting=True))
start = time.time()
stimulus = factory.random_rectangle_array(n_objects=800)
# assert isinstance(stimulus, nsn.RectangleArray)
# nsn.scale.log_size(stimulus, 1.2)
print(time.time()-start)

r = Rectangle(xy=(0, -40), size=(20, 20), attribute="blue")
r2 = Rectangle(xy=(-5, 50), size=(90, 160), attribute="red")

stim = pynsn.RectangleArray(target_area=pynsn.Dot(diameter=300))
stim.add([r2, r])

img = pil_image.create(stimulus, my_colours)
img.save("demo.png")

exit()

stimulus.mod_realign(keep_convex_hull=False, strict=False)
# stimulus.properties.fit_average_perimeter(130)

img = mpl_figure.create(stimulus, my_colours)
pp.savefig("demo2.png")
