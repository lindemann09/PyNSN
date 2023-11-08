__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math as _math
from multiprocessing import Pool as _Pool

import pygame as _pygame
from expyriment.misc import Clock as _Clock
from expyriment.stimuli import Canvas as _Canvas

from . import _base
from . import pil_image as _pil_image


class ExprimentNSNStimulus(_Canvas):

    def __init__(self, nsn_stimulus,
                 position=(0, 0),
                 antialiasing=True):

        _base.check_nsn_stimulus(nsn_stimulus)

        _Canvas.__init__(self, size=(0, 0), position=position)
        self.dot_array = nsn_stimulus
        self.antialiasing = antialiasing
        self._image = None

    @property
    def image(self):
        if self._image is None:
            self._image = _pil_image.create(nsn_stimulus=self.dot_array,
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
                 antialiasing=None,
                 make_pil_images_now=False,
                 multiprocessing=False):

        self.da_sequence = da_sequence
        self.stimuli = []
        self.position = position
        self.antialiasing = antialiasing

        for da in self.da_sequence.dot_arrays:
            stim = ExprimentNSNStimulus(nsn_stimulus=da, position=position,
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
