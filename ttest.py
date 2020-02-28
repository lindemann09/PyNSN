import pynsn
from pynsn import dot_array_sequence, Features, dot_array_archive
import json
import numpy as np

da_specs = pynsn.random_dot_array.Specs(
    target_area_radius=200,
    item_diameter_mean=5,
    item_diameter_range=(3, 8),
    item_diameter_std=2,
    item_colour=pynsn.Colour("skyblue"),
    minimum_gap=2)

seq = dot_array_sequence.create(specs=da_specs, match_feature=Features.ITEM_DIAMETER,
                           match_value=5,
                           min_max_numerosity=[10, 50],
                           round_decimals=0)


da = dot_array_archive.DotArrayArchive()
da.add(seq)
da.save("archive_demo.json", indent=2)

