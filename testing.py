from pynsn import factory, ImageColours, distr, RectangleArray, Colour, Rectangle, Coordinate2D, VisualFeature, Dot
from pynsn.image import pil, pyplot, svg
import numpy as np

a = Colour("red")
print(a.rgb)
print(a.rgb_alpha(1.0))

# define the visual features of the  dot array
da_specification2 = factory.DotArraySpecs(
    target_area_radius=200,
    diameter_distribution=distr.Beta(min_max=(10, 30), mu=15, sigma=2),
    minimum_gap=2)
da_specification = factory.RectangleArraySpecs(
    target_area_radius=200,
    width_distribution=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
    height_distribution=distr.Normal(min_max=(10, 40), mu=20, sigma=10),
    minimum_gap=2)
my_colours = ImageColours(target_area="#EEEEEE",
                          background=None,
                          item_colour="darkmagenta",
                          field_area_position="magenta",
                          field_area_outer = "blue",
                          #center_of_mass="red",
                          #center_of_outer_positions ="yellow",
                          )


stimulus = factory.random_array(da_specification, n_objects=15,
                                attributes=["blue", "green"])

# print(stimulus.split_array_by_attributes())
# print(stimulus._features.get_features_text())
# print(stimulus2.save("mystim.json", indent=2))

p = pil.create(stimulus, my_colours, antialiasing=True)
p.save("demo.png")

svg = svg.create(stimulus, my_colours, filename="demo.svg")
svg.save()



from matplotlib import pyplot as plt
f = pyplot.create(stimulus, my_colours)
plt.savefig("demo_pyplot.png", format="png")
#plt.show()