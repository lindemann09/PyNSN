"""
"""

from __future__ import unicode_literals, absolute_import, print_function
from builtins import zip, filter, range, super

from collections import namedtuple


RandomDotArrayImageAllParameter = namedtuple("DotPixmapParameter", # all paramter of da generator and pil image
                                "number max_array_radius dot_colour dot_diameter_mean " +
                                "dot_diameter_range dot_diameter_std minimum_gap colour_area " +
                                 "colour_convex_hull_positions colour_convex_hull_dots colour_center_of_mass "+
                                 "colour_center_of_outer_positions antialiasing colour_background")




ICON = RandomDotArrayImageAllParameter(number=11,
                           max_array_radius=200,
                           dot_colour="expyriment_orange",
                           dot_diameter_mean=35,
                           dot_diameter_range=[5, 80],
                           dot_diameter_std=20,
                           minimum_gap=2,
                           colour_area="#3e3e3e",
                           colour_convex_hull_positions=None,
                           colour_convex_hull_dots="expyriment_purple",
                           colour_center_of_mass=None,
                           colour_center_of_outer_positions=None,
                           antialiasing=True,
                           colour_background=None)
