from pynsn import factory, Colour, ImageColours, DotArray
from pynsn.image import svg


# define the visual features of the  dot array
da_specification = factory.DotArraySpecs(
    target_area_radius=200,
    item_diameter_mean=15,
    item_diameter_range=(10, 30),
    item_diameter_std=2,
    minimum_gap=2)

# generate on array with 100 dots
stimulus = factory.random_array(da_specification, 5)
stimulus2 = factory.random_array(da_specification, 4, attribute="skyblue",
                                 occupied_space=stimulus)
stimulus.join(stimulus2)



for x in stimulus.get():
    print(x)


#print(stimulus.split_array_by_attributes())
#print(stimulus._features.get_features_text())
#print(stimulus2.save("mystim.json", indent=2))
sv = svg.create(stimulus,  filename="demo.svg")
sv.save()



