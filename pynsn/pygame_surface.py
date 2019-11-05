__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import pygame

from .lib import pil_image

def create(dot_array,
           colours=pil_image.ImageColours(),
           antialiasing=True):

    if not isinstance(colours, pil_image.ImageColours):
        raise ValueError("Colours must be a ImageColours instance")

    img = pil_image.create(dot_array=dot_array,
                           colours=colours,
                           antialiasing=antialiasing)

    return pygame.image.fromstring(img.tobytes(), img.size, img.mode)
