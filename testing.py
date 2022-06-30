from pynsn import arrays, distr, random_array, adapt, scale, VisualFeature
from pynsn.image import svg_file, ImageColours

# define the visual features of the  dot array
ref = arrays.GenericObjectArray(target_area_radius=200)

size_dist = random_array.SizeDistribution(
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
                               size_distribution=size_dist,
                               n_objects=15,
                               attributes=["blue", "green"])

print(stimulus.features.log_size)
print(stimulus.features.log_size*1.2)

scale.visual_feature(stimulus, feature=VisualFeature.LOG_SIZE,
                                    factor=1.2)
print(stim_scaled.features.get(feature))


svg = svg_file.create(stim_matched, my_colours, filename="demo.svg")
svg.save()
