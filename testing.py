from pynsn import random_dot_array, Colour, ImageColours, DotArray
from pynsn.image import svg


# define the visual features of the  dot array
da_specification = random_dot_array.Specs(
    target_area_radius=200,
    item_diameter_mean=15,
    item_diameter_range=(10, 30),
    item_diameter_std=2,
    minimum_gap=2)

# generate on array with 100 dots
stimulus = random_dot_array.create(5, da_specification)
stimulus2 = random_dot_array.create(4, da_specification, attribute="skyblue",
                                    occupied_space=stimulus)
stimulus.join(stimulus2)



for x in stimulus.get():
    print(x)


#print(stimulus.split_array_by_attributes())
#print(stimulus._features.get_features_text())
#print(stimulus2.save("mystim.json", indent=2))
sv = svg.create(stimulus,  filename="demo.svg")
sv.save()



