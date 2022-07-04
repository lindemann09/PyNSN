import pynsn as nsn
from pynsn import distributions as distr
from pynsn.image import svg_file, pil_image, mpl_figure
from pynsn.visual_properties import fit

from pynsn import visual_properties as props
import numpy as np

# FIXME overlapping rects

euro_col = nsn.ImageColours(target_area="#003399", background="#003399")
euro_array = nsn.RectangleArray(target_area_radius=200, min_dist_between=1)
pict_file = nsn.PictureFile("eurostar.png")

for x in np.linspace(start=0, stop=2*np.pi, num=12+1):
    if x < 2*np.pi:
        pos = np.sin(x)*140, np.cos(x)*140
        pict = nsn.Rectangle(xy=pos, size=(40, 40), attribute=pict_file)
        euro_array.add(pict)

img = pil_image.create(euro_array, euro_col)
img.save("demo.png")


factory = nsn.NSNFactory(target_area_radius=200, min_dist_between=2)
factory.set_size_distribution(width=distr.Discrete((10, 20, 30)),
                              rectangle_proportion=1)
euro_array = factory.create_random_array(n_objects=27, attributes=pict_file)


# define the visual features of the  dot array
size_dist_dot = nsn.SizeDistribution(
    diameter=distr.Beta(min_max=(10, 30), mu=15, sigma=2)
)
size_dist_rect = nsn.SizeDistribution(
    width=distr.Normal(min_max=(10, 20), mu=15, sigma=2),
    rectangle_proportion=distr.Discrete([1, 2])
)

factory = nsn.NSNFactory(target_area_radius=150,
                         size_distribution=size_dist_rect)
stimulus = factory.create_random_array(n_objects=20,
                                       attributes=["blue", "green"])


####
stim_scaled = stimulus.copy()
#stim_scaled.center_array()
props.scale.log_size(stim_scaled, 1.5)


#stim_scaled.realign()
#props.scale.visual_property(stim_scaled, feature=feat, factor=1)
my_colours = nsn.ImageColours(target_area="#EEEEEE",
                          background=None,
                          opacity_object=0.9,
                          default_object_colour="darkmagenta",
                          field_area_positions="magenta",
                          field_area="blue",
                          center_of_field_area="red",
                          #center_of_mass="green"
                          )

img = pil_image.create(stimulus, my_colours)
img.save("demo_scaled.png")
svg = svg_file.create(stimulus, my_colours, filename="demo_scaled.svg")
svg.save()

mpl = mpl_figure.create(stimulus, my_colours)
mpl.savefig("demo_mpl.png")