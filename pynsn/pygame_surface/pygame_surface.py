
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import pygame

from ..pil_image import pil_image

# todo make PygameSurfaceGenerator

def create(dot_array,
           colour_target_area=None,
           colour_field_area=None,
           colour_field_area_outer=None,
           colour_center_of_mass = None,
           colour_center_of_outer_positions=None,
           antialiasing=True,
           colour_background=(0, 0, 0)):

    gen = pil_image.PILImageGenerator(colour_target_area=colour_target_area,
                                      colour_field_area=colour_field_area,
                                      colour_field_area_outer=colour_field_area_outer,
                                      colour_center_of_mass = colour_center_of_mass,
                                      colour_center_of_outer_positions=colour_center_of_outer_positions,
                                      antialiasing=antialiasing,
                                      colour_background=colour_background)

    img = gen.make(dot_array=dot_array)

    return pygame.image.fromstring(img.tobytes(), img.size, img.mode)
