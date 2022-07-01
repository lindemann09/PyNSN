from pynsn import arrays, random_array
from pynsn import distributions as distr
from pynsn.image import svg_file, ImageColours

from pynsn import visual_properties as props


# define the visual features of the  dot array
ref = arrays.GenericObjectArray(target_area_radius=200)

size_dist_dot = random_array.SizeDistribution(
    diameter=distr.Beta(min_max=(10, 30), mu=15, sigma=2)
)
size_dist_rect = random_array.SizeDistribution(
    width=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
    height=distr.Normal(min_max=(10, 40), mu=20, sigma=10)
)

my_colours = ImageColours(target_area="#EEEEEE",
                          background=None,
                          object_opacity=0.9,
                          default_object_colour="darkmagenta",
                          #field_area_positions="magenta",
                          field_area="blue",
                          #center_of_positions="red",
                          #center_of_mass="yellow"
                          )

stimulus = random_array.create(reference_array=ref,
                               size_distribution=size_dist_rect,
                               n_objects=5,
                               attributes=["blue", "green"])
props.scale.log_size(stimulus, 1.3)
svg = svg_file.create(stimulus, my_colours, filename="demo.svg")
svg.save()


dist = stimulus.distances_matrix(between_positions=False,
                                 overlap_is_zero=False)
print(dist<=0)
exit()

feat = props.flags.LOG_SIZE
stim_scaled = stimulus.copy()
props.scale.visual_feature(stim_scaled, feature=feat, factor=1)

svg = svg_file.create(stim_scaled, my_colours, filename="demo.svg")
svg.save()
props.fit.numerosity(stim_scaled, 10, keeping_field_area=True)
print(stim_scaled.features.field_area)
#stim_scaled.realign()
svg = svg_file.create(stim_scaled, my_colours, filename="demo_scaled.svg")
svg.save()
print(stim_scaled.features.field_area)
