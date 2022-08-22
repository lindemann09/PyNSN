"""_summary_"""

import time
from matplotlib import pyplot as pp
import numpy as np
from numpy.typing import NDArray, ArrayLike

import pynsn
import pynsn as nsn
from pynsn import distributions as distr
from pynsn._shapes.rectangle import Rectangle
from pynsn.image import svg_file, pil_image, mpl_figure
from pynsn._lib import np_rectangles, np_dots, np_coordinates, np_tools

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


a = np.array((0, 0))
b = np.array(((-10, -10), (-100, 0), (103, 303)))
s = np.array([100, 100])

polar = np_coordinates.cartesian2polar(
    np.asarray(b) - np.asarray(a))  # type: ignore

x = np_rectangles._center_edge_distance(polar[:, 1], s)
# x = required_displacement(a, s, b, s)
print(x)

exit()

factory = nsn.NSNFactory(min_distance_between_objects=2,
                         # target_area=pynsn.Rectangle(size=(200, 400)),
                         target_area=pynsn.Dot(diameter=400),
                         min_distance_area_boarder=2)

factory.set_appearance_dots(diameter=(20, 30, 40),
                            attributes=distr.Levels(["blue", "green"],
                                                    exact_weighting=True))

factory.set_appearance_rectangles(width=(20, 30, 40), proportion=1,
                                  attributes=distr.Levels(["blue", "green"],
                                                          exact_weighting=True))


# start = time.time()
stimulus = factory.random_rectangle_array(n_objects=20)
# assert isinstance(stimulus, nsn.RectangleArray)
# nsn.scale.log_size(stimulus, 1.2)
# print(time.time()-start)

img = pil_image.create(stimulus, my_colours)
img.save("demo.png")

# print(stimulus)


# stimulus.mod_realign(keep_convex_hull=False, strict=False)

stimulus.properties.fit_field_area(99000)
# print(stimulus)
# stimulus.mod_squeeze_to_area()  # FIXME squeed too much, because incorrect distance calulations

idx = stimulus.get_outlier()

img = mpl_figure.create(stimulus, my_colours)
pp.savefig("demo2.png")
