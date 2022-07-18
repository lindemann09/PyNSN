import pynsn as nsn
from pynsn import distributions as distr
from pynsn.image import svg_file, pil_image, mpl_figure
from pynsn.visual_properties import fit
from time import monotonic
from pynsn import visual_properties as props
import numpy as np

from random import randint
seed = randint(0, 1000)
seed = 364
nsn.init_random_generator(seed)

# FIXME overlapping rects
my_colours = nsn.ImageColours(#target_area="#EEEEEE",
                              background=None,
                              opacity_object=0.9,
                              default_object_colour="darkmagenta",
                              #field_area_positions="magenta",
                              field_area="gray",
                              #center_of_field_area="red",
                              # center_of_mass="green"
                              )


factory = nsn.NSNFactory(target_area_radius=200, min_dist_between=2)
factory.set_appearance_dot(diameter=(40, 10, 30), attributes=distr.Levels(["blue", "green"],
                                        exact_weighting=True) )

stimulus = factory.create_random_array(n_objects=10)
props.scale.log_size(stimulus, 1.4)

img = pil_image.create(stimulus, my_colours)
img.save("demo.png")

####
print("start")
t = monotonic()
stimulus.remove_overlaps(keep_field_area=True, strict=False)
# props.scale.visual_property(stim_scaled, feature=feat, factor=1)
print("time {}".format((monotonic()-t)))

img = pil_image.create(stimulus, my_colours)
img.save("demo_scaled.png")
