from pynsn import random_dot_array, Colour, ImageColours, DotArray
from pynsn.image import svg

c = Colour(None)

# define the visual features of the  dot array
da_specification = random_dot_array.Specs(
    target_area_radius=200,
    item_diameter_mean=15,
    item_diameter_range=(10, 30),
    item_diameter_std=2,
    minimum_gap=2)

# generate on array with 100 dots
stimulus = random_dot_array.create(100, da_specification, attributes="yellow")
#stimulus.set_attributes("red")
print(stimulus.features.get_features_text())

#stimulus.save("mystim.json")

da = DotArray(0,0)
da.load("mystim.json")
sv = svg.create(da, colours=ImageColours(default_dot_colour="blue"), filename="demo.svg")
sv.save()


