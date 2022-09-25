"""_summary_"""

from pyarrow import feather
import pyarrow
import pandas
import timeit
from matplotlib import pyplot as pp
import numpy as np
from numpy.typing import NDArray, ArrayLike

import pynsn
import pynsn as nsn
from pynsn import distributions as distr
from pynsn.spatial_relations import _spatial_relations
from pynsn._shapes.rectangle import Rectangle
from pynsn.image import svg_file, pil_image, mpl_figure
from pynsn._lib import geometry, np_tools

seed = 921
nsn.init_random_generator(seed)

data = pynsn.NSNStimulus(object_array=nsn.DotArray(),
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


def timetest():
    factory.random_rectangle_array(n_objects=20)


#print("create:", timeit.timeit(timetest, number=100)*10)


stimulus = factory.random_rectangle_array(n_objects=20)


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

a = pynsn.Picture("examples/eurostar.png")
# stimulus.objects.set_attributes(a.attribute)
img = pil_image.create(stimulus, my_colours)
img.save("demo.png")

# print(stimulus)


# stimulus.mod_realign(keep_convex_hull=False, strict=False)


stimulus.properties.fit_field_area(99000)
print(stimulus)

idx = stimulus.get_outlier()

# stimulus.mod_squeeze_to_area() # FIXME squeed does not work

img = pil_image.create(stimulus, my_colours)


def timetest2():
    pil_image.create(stimulus, my_colours)
#print("picture:", timeit.timeit(timetest2, number=100)*10)


img.save("demo2.png")
