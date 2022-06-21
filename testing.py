from pynsn import factory, ImageColours, distr
from pynsn.image import pil

import numpy as np




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
