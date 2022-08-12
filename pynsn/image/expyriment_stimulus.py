__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math as _math
from multiprocessing import Pool as _Pool

import pygame as _pygame
from expyriment.misc import Clock as _Clock
from expyriment.stimuli import Canvas as _Canvas

from . import _array_draw, _colour
from . import pil_image as _pil_image


class ExprimentDotArray(_Canvas):

    def __init__(self, object_array,
                 colours=None,
                 position=(0, 0),
                 antialiasing=True):

        _array_draw._check_object_array(object_array)
        if colours is None:
            colours = _colour.ImageColours()
        if not isinstance(colours, _colour.ImageColours):
            raise TypeError("Colours must be of type image.ImageColours")

        _Canvas.__init__(self, size=(0, 0), position=position)
        self.dot_array = object_array
        self.colours = colours
        self.antialiasing = antialiasing
        self._image = None

    @property
    def image(self):
        if self._image is None:
            self._image = _pil_image.create(object_array=self.dot_array,
                                            colours=self.colours,
                                            antialiasing=self.antialiasing)  # TODO gabor filter

        return self._image

    def _create_surface(self):
        self._size = self.image.size
        return _pygame.image.frombuffer(self.image.tobytes(),
                                        self.image.size,
                                        self.image.mode)


class ExpyrimentDASequence(object):

    def __init__(self, da_sequence,
                 # pil_image_generator TODO better using generator
                 position=(0, 0),
                 colours=_colour.ImageColours(),
                 antialiasing=None,
                 make_pil_images_now=False,
                 multiprocessing=False):

        if not isinstance(colours, _colour.ImageColours):
            raise TypeError("Colours must be a ImageColours instance")

        self.da_sequence = da_sequence
        self.stimuli = []
        self.position = position
        self.antialiasing = antialiasing

        for da in self.da_sequence.dot_arrays:
            stim = ExprimentDotArray(object_array=da, position=position,
                                     colours=colours,
                                     antialiasing=antialiasing)
            self.stimuli.append(stim)

        if make_pil_images_now:

            if not multiprocessing:
                list(map(lambda x: x._create_pil_image(), self.stimuli))
                self._make_image_process = None
            else:
                p = _Pool()

                for c, pil_im in enumerate(p.imap(ExpyrimentDASequence._make_stimuli_map_helper, self.stimuli)):
                    self.stimuli[c]._image = pil_im
                p.close()
                p.join()

    def get_stimulus_numerosity(self, number_of_dots):
        """returns image with a particular numerosity"""
        try:
            return self.stimuli[self.da_sequence.numerosity_idx[number_of_dots]]
        except IndexError:
            return None

    @property
    def is_preloaded(self):
        for x in reversed(self.stimuli):
            if not x.is_preloaded:
                return False
        return True

    def preload(self, until_percent=100, time=None, do_not_return_earlier=False):
        """
        preloaded all _lib stimuli

        Note: this will take a while!

        preload certain percent or or a time.

        """
        if until_percent > 0 and until_percent < 100:
            last = int(_math.floor(until_percent * len(self.stimuli) / 100.0))
        elif until_percent == 0:
            last = 0
        else:
            last = len(self.stimuli)
        cl = _Clock()

        try:
            for x in self.stimuli[:last]:
                if not x.is_preloaded and (time is None or cl.time < time):
                    x.preload()
            rtn = True
        except:
            rtn = False

        if do_not_return_earlier and time is not None:
            cl.wait(time - cl.time)
        return rtn

    def unload(self):
        """
        returns array of preloaded dot_array_sequence
        """

        try:
            list(map(lambda x: x.unload(), self.stimuli))
            return True
        except:
            return False

    @staticmethod
    def _make_stimuli_map_helper(x):
        return x._create_pil_image()
