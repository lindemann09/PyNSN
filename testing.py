from pynsn import random_dot_array, Colour
from pynsn.image import svg

# define the visual features of the  dot array
da_specification = random_dot_array.Specs(
    target_area_radius=200,
    item_diameter_mean=15,
    item_diameter_range=(10, 30),
    item_diameter_std=2,
    item_colour=Colour("skyblue"),
    minimum_gap=2)

# generate on array with 100 dots
rda = random_dot_array.create(100, da_specification)
print(rda.features.get_features_text())

sv = svg.create(rda, filename="demo.svg")
sv.save()
print(rda.json(include_hash=True))

