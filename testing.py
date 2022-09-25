"""_summary_"""

from pyarrow import feather
import pyarrow
import pandas
import time
from matplotlib import pyplot as pp
import numpy as np
from numpy.typing import NDArray, ArrayLike

import pynsn
import pynsn as nsn
from pynsn import distributions as distr
from pynsn._shapes import spatial_relations
from pynsn._shapes.rectangle import Rectangle
from pynsn.image import svg_file, pil_image, mpl_figure
from pynsn._lib import geometry, np_tools

seed = 921
nsn.init_random_generator(seed)

dota = pynsn.NSNStimulus(object_list=nsn.DotList(),
                         target_area_shape=pynsn.Dot(diameter=400))

factory = nsn.NSNFactory(min_dist_between_objects=2,
                         # target_area=pynsn.Rectangle(size=(200, 400)),
                         target_area=pynsn.Dot(diameter=400),
                         min_dist_area_edge=0)

factory.set_appearance_dots(diameter=(20, 30, 40),
                            attributes=distr.Levels(["blue", "green"],
                                                    exact_weighting=True))

factory.set_appearance_rectangles(width=(21, 30, 40), proportion=1,
                                  attributes=distr.Levels(["blue", "green"],
                                                          exact_weighting=True))


# start = time.time()
stimulus = factory.random_rectangle_array(n_objects=20)
# assert isinstance(stimulus, nsn.RectangleArray)
# nsn.scale.log_size(stimulus, 1.2)
# print(time.time()-start)


# make image

my_colours = pil_image.ImageColours(  # target_area="#EEEEEE",
    background="",
    opacity_object=0.9,
    default_object_colour="darkmagenta",
    # field_area_positions="magenta",
    # field_area="gray",
    # center_of_positions="red",
    # center_of_mass="magenta"
)


img = pil_image.create(stimulus, my_colours)
img.save("demo.png")

# print(stimulus)


# stimulus.mod_realign(keep_convex_hull=False, strict=False)

stimulus.properties.fit_field_area(99000)
# print(stimulus)

idx = stimulus.get_outlier()
print(stimulus)
# FIXME squeed too much, because incorrect distance calulations
# stimulus.mod_squeeze_to_area()

img = mpl_figure.create(stimulus, my_colours)
pp.savefig("demo2.png")
