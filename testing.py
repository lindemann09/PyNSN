import pynsn
from pynsn import VisualFeatures
from pynsn.dot_array import random_dot_array
from pynsn.sequence import dot_array_archive, dot_array_sequence

da_specs = pynsn.dot_array.random_dot_array.Specs(
    target_area_radius=200,
    item_diameter_mean=5,
    item_diameter_range=(3, 8),
    item_diameter_std=2,
    item_colour=pynsn.Colour("skyblue"),
    minimum_gap=2)


rda = random_dot_array.create(5, da_specs)
#print(rda.convex_hull.scipy_convex_hull)
print(rda.features.field_area)
print(rda.features.mean_item_diameter)
print(rda.features.da.xy)

rda.delete(2)
rda.delete(3)

#print(rda.convex_hull.xy)
print(rda.features.da.xy)
print(rda.features.field_area)
print(rda.features.mean_item_diameter)

exit()









seq = dot_array_sequence.create(specs=da_specs, match_feature=VisualFeatures.ITEM_DIAMETER,
                                match_value=5,
                                min_max_numerosity=[10, 50],
                                round_decimals=0)



da = dot_array_archive.DotArrayArchive()
da.add(seq)
da.save("archive_demo.json", indent=2)

