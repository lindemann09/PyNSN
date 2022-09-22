__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Optional, Union

import pygame as _pygame

from .. import _arrays
from . import _array_draw
from . import pil_image as _pil_image
from ._image_colours import ImageColours


def create(object_array: _arrays.ObjectArrayType,
           colours: Optional[ImageColours] = None,
           antialiasing: Union[bool, int] = True) -> _pygame.Surface:

    _array_draw._check_object_array(object_array)
    if colours is None:
        colours = ImageColours()
    if not isinstance(colours, ImageColours):
        raise TypeError("Colours must be of type image.ImageColours")

    img = _pil_image.create(object_array=object_array,
                            colours=colours,
                            antialiasing=antialiasing)

    return _pygame.image.fromstring(img.tobytes(), img.size, img.mode)
