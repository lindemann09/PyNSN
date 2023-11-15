__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Union

import pygame as _pygame

from .. import _stimulus
from . import pil_image as _pil_image


def create(nsn_stimulus: _stimulus.NSNStimulus,
           antialiasing: Union[bool, int] = True) -> _pygame.Surface:

    img = _pil_image.create(nsn_stimulus=nsn_stimulus,
                            antialiasing=antialiasing)

    return _pygame.image.fromstring(img.tobytes(), img.size,
                                    img.mode)  # type: ignore
