import pynsn
import os
from pynsn import VisualFeatures, dot_array, ImageColours
from pynsn.dot_array import random_dot_array
from pynsn.sequence import dot_array_archive, dot_array_sequence


da_specs = random_dot_array.Specs(
    target_area_radius=200,
    item_diameter_mean=15,
    item_diameter_range=(10, 30),
    item_diameter_std=2,
    item_colour=pynsn.Colour("skyblue"),
    minimum_gap=2)

rda = pynsn.random_dot_array.create(100, da_specs)
rda.round(decimals=1)
rda.save("demp.json", indent=2, include_hash=False)

