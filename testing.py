from pynsn import factory, ImageColours, distr, RectangleArray, Rectangle
from pynsn.image import pil

import numpy as np

d = distr.Discrete(population=[1,2,3], weights=[10, 43, 10])
b = d.sample(10)
print(d.as_dict())


exit()
a = distr.Normal(min_max=(10, 100), mu=55, sigma=20)

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
my_colours = ImageColours(target_area="#EEEEEE", background="gray",
                          item_colour="darkmagenta")  # show named colours see Colour.NAMED_COLOURS


ra = RectangleArray(target_array_radius = 200,
                          minimum_gap = 0)
# generate on array with 100 dots
stimulus = factory.random_array(da_specification, n_objects=5,
                                attributes=["blue", "green"])

print(stimulus.features.as_dict())

# print(stimulus.split_array_by_attributes())
# print(stimulus._features.get_features_text())
# print(stimulus2.save("mystim.json", indent=2))

pil.create(stimulus, my_colours).save("demo.png")
exit()
