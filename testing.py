import pynsn as nsn
from pynsn import distributions as distr
from pynsn.image import svg_file, pil_image, mpl_figure
from pynsn.visual_properties import fit
from time import monotonic
from pynsn import visual_properties as props
import numpy as np

from random import randint
seed = randint(0, 1000)
seed = 993
nsn.init_random_generator(seed)

# FIXME overlapping rects
my_colours = nsn.ImageColours(#target_area="#EEEEEE",
                              background=None,
                              opacity_object=0.9,
                              default_object_colour="darkmagenta",
                              #field_area_positions="magenta",
                              field_area="gray",
                              #center_of_positions="red",
                              #center_of_mass="magenta"
                              )


factory = nsn.NSNFactory(target_area_radius=200, min_dist_between=2)
factory.set_appearance_dot(diameter=(40, 10, 30), attributes=distr.Levels(["blue", "green"],
                                        exact_weighting=True) )

stimulus = factory.create_random_array(n_objects=20)
assert isinstance(stimulus, nsn.DotArray)
props.scale.log_size(stimulus, 1.2)

img = pil_image.create(stimulus, my_colours)
img.save("demo.png")


stimulus.mod_remove_overlaps(keep_field_area=True, strict=False)
#stimulus.mod_center_field_area()
#stimulus.mod_squeeze_to_area()

img = pil_image.create(stimulus, my_colours)
img.save("demo_scaled.png")
