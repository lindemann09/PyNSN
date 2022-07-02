__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import pygame as _pygame

from . import _colour
from . import pil_image as _pil_image
from ..arrays import _check_object_array

from ._colour import ImageColours # make available

def create(object_array,
           colours=None,
           antialiasing=True):
    _check_object_array(object_array)
    if colours is None:
        colours = _colour.ImageColours()
    if not isinstance(colours, _colour.ImageColours):
        raise TypeError("Colours must be of type image.ImageColours")

    img = _pil_image.create(object_array=object_array,
                            colours=colours,
                            antialiasing=antialiasing)

    return _pygame.image.fromstring(img.tobytes(), img.size, img.mode)
