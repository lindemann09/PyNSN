from pynsn import arrays, distr, random_array
from pynsn.image import pil_image, mpl_figure, svg_file, ImageColours
import numpy as np


# define the visual features of the  dot array
ref = arrays.GenericObjectArray(target_area_radius=200)

size_dist = random_array.SizeDistribution(
    diameter=distr.Beta(min_max=(10, 30), mu=15, sigma=2)
)
size_dist2 = random_array.SizeDistribution(
    width=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
    height=distr.Normal(min_max=(10, 40), mu=20, sigma=10)
)

my_colours = ImageColours(target_area="#EEEEEE",
                          background=None,
                          object_opacity=0.9,
                          default_object_colour="darkmagenta",
                          field_area_positions="magenta",
                          field_area="blue",
                          center_of_positions="red",
                          center_of_mass="yellow",
                          info_shapes_opacity=0.5/3
                          )



stimulus = random_array.create(reference_array=ref,
                               size_distribution=size_dist2,
                               n_objects=15,
                               attributes=["blue", "green"])

p = pil_image.create(stimulus, my_colours, antialiasing=True)
p.save("demo.png")

svg = svg_file.create(stimulus, my_colours, filename="demo.svg")
svg.save()

from matplotlib import pyplot as plt
f = mpl_figure.create(stimulus, my_colours)
plt.savefig("demo_pyplot.png", format="png")
#plt.show()