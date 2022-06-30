from pynsn import arrays, distr, random_array, adapt, scale, VisualFeatureTypes
from pynsn.image import svg_file, ImageColours

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
                               size_distribution=size_dist_dot,
                               n_objects=15,
                               attributes=["blue", "green"])


feat = VisualFeatureTypes.LOG_SIZE
stim_scaled = stimulus.copy()
scale.visual_feature(stim_scaled, feature=feat, factor=1)

svg = svg_file.create(stim_scaled, my_colours, filename="demo.svg")
svg.save()
adapt.numerosity(stim_scaled, 10, keeping_field_area=True)
print(stim_scaled.features.field_area)
#stim_scaled.realign()
svg = svg_file.create(stim_scaled, my_colours, filename="demo_scaled.svg")
svg.save()
print(stim_scaled.features.field_area)
