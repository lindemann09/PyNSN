"""_summary_"""

import numpy as np
from numpy.typing import NDArray, ArrayLike

import pynsn
import pynsn as nsn
from pynsn import distributions as distr
from pynsn.image import svg_file, pil_image, mpl_figure

from random import randint
import numpy as np


seed = 991
nsn.init_random_generator(seed)

my_colours = nsn.ImageColours(  # target_area="#EEEEEE",
    background=None,
    opacity_object=0.9,
    default_object_colour="darkmagenta",
    # field_area_positions="magenta",
    field_area="gray",
    # center_of_positions="red",
    # center_of_mass="magenta"
)


r = pynsn.Rectangle(xy=(0, 0), size=(80, 80),
                    attribute=nsn.PictureFile("mypict.png"))


factory = nsn.NSNFactory(min_distance_between_objects=2,
                         target_area=pynsn.Rectangle(size=(200, 400)),
                         # starget_area=pynsn.Dot(diameter=400),
                         min_distance_area_boarder=2)

factory.set_appearance_dot(diameter=(10, 10, 10),
                           attributes=distr.Levels(["blue", "green"],
                                                   exact_weighting=True))

# factory.set_appearance_rectangle(width=(10, 10, 10), proportion=0.5,
#                                 attributes=distr.Levels(["blue", "green"],
#                                                         exact_weighting=True))

stimulus = factory.create_random_array(n_objects=190)
# assert isinstance(stimulus, nsn.RectangleArray)
# nsn.scale.log_size(stimulus, 1.2)

img = pil_image.create(stimulus, my_colours)
img.save("demo.png")

exit()
stimulus.mod_realign(keep_convex_hull=False, strict=False)

stimulus.properties.fit_average_perimeter(130)
dist = stimulus.get_distances_matrix()

img = pil_image.create(stimulus, my_colours)
img.save("demo_scaled.png")
