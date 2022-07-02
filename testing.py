from pynsn import arrays, random_array
from pynsn import distributions as distr
from pynsn.image import svg_file, ImageColours
from pynsn import shapes

# FIXME Own exception
# FIXME __iter__ and iter_objects


from pynsn import visual_properties as props
import numpy as np

# define the visual features of the  dot array
ref = arrays.DotArray(target_area_radius=200)

size_dist_dot = random_array.SizeDistribution(
    dot_diameter=distr.Beta(min_max=(10, 30), mu=15, sigma=2)
)
size_dist_rect = random_array.SizeDistribution(
    rectangle_width=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
    #rectangle_height=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
    rectangle_proportion=distr.Discrete([1, 2])
)


my_colours = ImageColours(target_area="#EEEEEE",
                          background=None,
                          object_opacity=0.9,
                          default_object_colour="darkmagenta",
                          field_area_positions="magenta",
                          field_area="blue",
                          center_of_field_area="red",
                          #center_of_mass="green"
                          )

stimulus = random_array.create(reference_array=ref,
                               size_distribution=size_dist_dot,
                               n_objects=7,
                               attributes=["blue", "green"])
props.scale.log_size(stimulus, 1.3)
svg = svg_file.create(stimulus, my_colours, filename="demo.svg")
svg.save()



a = stimulus.iter_objects()
b = list(a)

exit()

####
stim_scaled = stimulus.copy()
#stim_scaled.center_array()


#stim_scaled.realign()
#props.scale.visual_property(stim_scaled, feature=feat, factor=1)
svg = svg_file.create(stim_scaled, my_colours, filename="demo_scaled.svg")
svg.save()
