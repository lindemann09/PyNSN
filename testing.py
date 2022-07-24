import pynsn as nsn
from pynsn import distributions as distr
from pynsn.image import svg_file, pil_image, mpl_figure
import numpy as np
from random import randint
seed = randint(0, 1000)
seed = 991
nsn.init_random_generator(seed)

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

#factory.set_appearance_rectangle(width=(40, 10, 30), proportion=0.5,
#                                        attributes=distr.Levels(["blue", "green"],
#                                        exact_weighting=True))

stimulus = factory.create_random_array(n_objects=20)
stimulus.copy(indices=np.asarray([3,5,6]))


#assert isinstance(stimulus, nsn.RectangleArray)
#nsn.scale.log_size(stimulus, 1.2)

img = pil_image.create(stimulus, my_colours)
img.save("demo.png")
stimulus.mod_realign(keep_convex_hull=False, strict=False)

stimulus.properties.fit_average_perimeter(130)
dist = stimulus.get_distances_matrix()

img = pil_image.create(stimulus, my_colours)
img.save("demo_scaled.png")

