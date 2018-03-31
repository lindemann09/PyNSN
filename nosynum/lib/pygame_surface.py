from __future__ import absolute_import, print_function, division

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import pygame
from . import pil_image

def create(dot_array,
           colour_area=None,
           colour_convex_hull_positions=None,
           colour_convex_hull_dots=None,
           colour_center_of_mass = None,
           colour_center_of_outer_positions=None,
           antialiasing=None,
           colour_background=(0, 0, 0)):

    img = pil_image.create(dot_array=dot_array,
                    colour_area=colour_area,
                    colour_convex_hull_positions=colour_convex_hull_positions,
                    colour_convex_hull_dots=colour_convex_hull_dots,
                    colour_center_of_mass = colour_center_of_mass,
                    colour_center_of_outer_positions=colour_center_of_outer_positions,
                    antialiasing=antialiasing,
                    colour_background=colour_background)

    return pygame.image.fromstring(img.tobytes(), img.size, img.mode)
