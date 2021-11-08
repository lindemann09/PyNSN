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


seq = dot_array_sequence.create(specs=da_specs, match_feature=VisualFeatures.ITEM_DIAMETER,
                                match_value=5,
                                min_max_numerosity=[10, 50],
                                round_decimals=0)



da = dot_array_archive.DotArrayArchive()
da.add(seq)
da.save("archive_demo.json", indent=2)

