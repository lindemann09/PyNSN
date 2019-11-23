

import pynsn
from pynsn.expyriment_stimulus import ExprimentDotArray
import json

da_specs = pynsn.random_dot_array.Specs(
    target_area_radius=200,
    item_diameter_mean=5,
    item_diameter_range=(3, 8),
    item_diameter_std=2,
    item_colour=pynsn.Colour("skyblue"),
    minimum_gap=2)


stim = pynsn.random_dot_array.create(n_dots=32, specs=da_specs)
stim.round(decimals=0)
stim.save("/tmp/test.json")

stim2 = pynsn.DotArray(target_array_radius=100, minimum_gap=2)
stim2.load("/tmp/test.json")
print(stim2.xy)
print(stim2.get_attributes())
