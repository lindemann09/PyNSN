__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import pygame as _pygame

from . import _colour
from . import pil as _pil_image

def create(dot_array,
           colours=_colour.ImageColours(),
           antialiasing=True):

    if not isinstance(colours, _colour.ImageColours):
        raise ValueError("Colours must be a ImageColours instance")

    img = _pil_image.create(dot_array=dot_array,
                           colours=colours,
                           antialiasing=antialiasing)

    return _pygame.image.fromstring(img.tobytes(), img.size, img.mode)
