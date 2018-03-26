from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import pygame
from PIL.Image import isImageType
from expyriment.stimuli import Canvas

class ExprimentPILImage(Canvas):

    def __init__(self, pil_image, position=(0,0)):

        if not isImageType(pil_image):
            raise RuntimeError("pil_image (type: {}) is not a PIL image".format(type(pil_image)))

        Canvas.__init__(self, size=pil_image.size, position=position)
        self._image = pil_image

    def _create_surface(self):
        return pygame.image.fromstring(self._image.tobytes(),
                                       self.size,
                                       self._image.mode)