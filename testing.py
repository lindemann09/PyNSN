<<<<<<< HEAD
from pynsn import factory, ImageColours, distr
from pynsn.image import pil

import numpy as np

from pynsn.data_types import *


def scale(scalar: float, vector: NumVector, x : NumArray) -> NumVector:
    return [scalar * num for num in vector]


# typechecks; a list of floats qualifies as a Vector.
x = np.array([[98, 98], [89, 78]])
mat = [[1, 3, 2], [5, 6, 4]]
new_vector = scale(2.0, [1.0, "df", -4.2, 5.4], x=x)

a = distr.Normal(min_max=(10, 100), mu=55, sigma=20)

# define the visual features of the  dot array
da_specification2 = factory.DotArraySpecs(
    target_area_radius=200,
    diameter_distribution=distr.Beta(min_max=(10, 30), mu=15, sigma=2),
    minimum_gap=2)
da_specification = factory.RectangleArraySpecs(
    target_area_radius=200,
    width_distribution=distr.Normal(min_max=(10, 80), mu=30, sigma=10),
    height_distribution=distr.Normal(min_max=(10, 80), mu=30, sigma=10),
    minimum_gap=2)
my_colours = ImageColours(target_area="#EEEEEE", background="gray",
                          item_colour="darkmagenta")  # show named colours see Colour.NAMED_COLOURS

# generate on array with 100 dots
stimulus = factory.random_array(da_specification, 4)
stimulus.round(2)
# r = Rectangle(xy=(-85.76060604630113, -77.51417327204257), size=(20, 40), attribute=None)
# print(r)


for x in stimulus.get():
    print(x)

# print(stimulus.split_array_by_attributes())
# print(stimulus._features.get_features_text())
# print(stimulus2.save("mystim.json", indent=2))

pil.create(stimulus, my_colours).save("demo.png")
exit()
=======
import pynsn as nsn
from pynsn import distributions as distr
from pynsn.image import svg_file, pil_image, mpl_figure
from pynsn.visual_properties import fit
from time import monotonic
from pynsn import visual_properties as props
import numpy as np
from pynsn._lib import misc
from random import randint
seed = randint(0, 1000)
seed = 997
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

factory.set_appearance_rectangle(width=(40, 10, 30), proportion=0.5,
                                        attributes=distr.Levels(["blue", "green"],
                                        exact_weighting=True))

stimulus = factory.create_random_array(n_objects=20)
assert isinstance(stimulus, nsn.RectangleArray)
props.scale.log_size(stimulus, 1.2)

img = pil_image.create(stimulus, my_colours)
img.save("demo.png")
stimulus.mod_realign(keep_convex_hull=False, strict=False)
print(stimulus.save("demp.json"))

dist = stimulus.get_distances_matrix()

img = pil_image.create(stimulus, my_colours)
img.save("demo_scaled.png")

>>>>>>> devel
