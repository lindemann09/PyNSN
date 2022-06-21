__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import pygame as _pygame

from . import _colour
from . import pil as _pil_image
from .._lib.arrays import GenericObjectArray as _GenericObjectArray

def create(object_array,
           colours,
           antialiasing=True):

    assert isinstance(object_array, _GenericObjectArray)

    if not isinstance(colours, _colour.ImageColours):
        raise TypeError("Colours must be a ImageColours instance.")

    img = _pil_image.create(object_array=object_array,
                            colours=colours,
                            antialiasing=antialiasing)

    return _pygame.image.fromstring(img.tobytes(), img.size, img.mode)
