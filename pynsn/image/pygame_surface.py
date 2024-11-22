__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import pygame as _pygame

from . import pil_image as _pil_image
from .. import NSNStimulus as _NSNStimulus


def create(nsn_stimulus: _NSNStimulus,
           antialiasing: bool | int = True) -> _pygame.Surface:

    img = _pil_image.create(nsn_stimulus=nsn_stimulus,
                            antialiasing=antialiasing)

    return _pygame.image.fromstring(img.tobytes(), img.size,
                                    img.mode)  # type: ignore
