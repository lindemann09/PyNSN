import pynsn as nsn
from pynsn import random_array, RectangleArray, DotArray, Rectangle, Dot
from pynsn import distributions as distr
from pynsn.image import svg_file, ImageColours, pil_image, defaults
from pynsn.visual_properties import fit
from pynsn.distributions import Discrete
from pynsn._lib.shapes import picture_attr

from pynsn import visual_properties as props
import numpy as np

# FIXME test Attribute Array (properties)

euro_col = ImageColours(target_area="#003399", background="#003399")
euro_array = RectangleArray(target_area_radius=200, min_dist_between=1)

for x in np.linspace(start=0, stop=2*np.pi, num=12+1):
    if x < 2*np.pi:
        pos = np.sin(x)*140, np.cos(x)*140
        pict = Rectangle(xy=pos, size=(40, 40), picture="eurostar.png")
        euro_array.add(pict)

img = pil_image.create(euro_array, euro_col)
img.save("demo.png")


euro_array.shuffle_all_positions()
img = pil_image.create(euro_array, euro_col)
#img.save("demo.png")


pict_file= ""
size_dist = random_array.SizeDistribution(width=distr.Discrete((10, 20, 30)),
                                          rectangle_proportion=1)
euro_array = random_array.create(euro_array, size_dist, n_objects=27,
                                 attributes=picture_attr(pict_file))



s = (40, 40)
objs =  [Rectangle(xy=(-65, 50), size=s, picture=pict_file),
                Rectangle(xy=(25, 45), size=s, attribute="black"),
                Rectangle(xy=(60, -50), size=s, picture=pict_file),
                Rectangle(xy=(-80, -20), size=s, picture=pict_file),
                Rectangle(xy=(-0, 0), size=s, picture=pict_file)]
# euro_array = RectangleArray(target_area_radius=150, min_dist_between=3)
# euro_array.add(objs)
# euro_col = ImageColours(target_area="#003399")
# img = pil_image.create(euro_array, euro_col)
# img.save("demo.png")

exit()

# define the visual features of the  dot array


size_dist_dot = random_array.SizeDistribution(
    diameter=distr.Beta(min_max=(10, 30), mu=15, sigma=2)
)

size_dist_rect = random_array.SizeDistribution(
    width=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
    #rectangle_height=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
    rectangle_proportion=distr.Discrete([1, 2])
)

ref = RectangleArray(target_area_radius=200)

my_colours = ImageColours(target_area="#EEEEEE",
                          background=None,
                          opacity_object=0.9,
                          default_object_colour="darkmagenta",
                          field_area_positions="magenta",
                          field_area="blue",
                          center_of_field_area="red",
                          #center_of_mass="green"
                          )

stimulus = random_array.create(n_objects=7,
                               reference_array=ref,
                               size_distribution=size_dist_rect,
                               attributes=["blue", "green"])

props.scale.log_size(stimulus, 1.15)
#svg = svg_file.create(stimulus, my_colours, filename="demo.svg")
#svg.save()

incr = random_array.create_incremental(n_objects=10,
                                              reference_array=ref,
                                              size_distribution=size_dist_rect,
                                              attributes=["blue", "green"])

for c, arr in enumerate(incr):
    svg = svg_file.create(arr, my_colours, filename="demo{}.svg".format(c))
    svg.save()
exit()

print(stimulus.properties.convex_hull._convex_hull.simplices)
print(stimulus.properties.convex_hull._convex_hull.vertices)

####
stim_scaled = stimulus.copy()
#stim_scaled.center_array()


#stim_scaled.realign()
#props.scale.visual_property(stim_scaled, feature=feat, factor=1)
svg = svg_file.create(stim_scaled, my_colours, filename="demo_scaled.svg")
svg.save()
