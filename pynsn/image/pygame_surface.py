__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Optional, Union

import pygame as _pygame

from .. import _stimulus
from . import pil_image as _pil_image
from ._array_draw import check_nsn_stimulus
from ._image_colours import ImageColours

# FIXME check pygame2 compatibility


def create(nsn_stimulus: _stimulus.NSNStimulus,
           colours: Optional[ImageColours] = None,
           antialiasing: Union[bool, int] = True) -> _pygame.Surface:

    check_nsn_stimulus(nsn_stimulus)
    if colours is None:
        colours = ImageColours()
    if not isinstance(colours, ImageColours):
        raise TypeError("Colours must be of type image.ImageColours")

    img = _pil_image.create(nsn_stimulus=nsn_stimulus,
                            colours=colours,
                            antialiasing=antialiasing)

    return _pygame.image.fromstring(img.tobytes(), img.size, img.mode)
