from __future__ import absolute_import
from builtins import zip, filter, range, super, map

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math
import pygame
from multiprocessing import Pool

from expyriment.misc import Clock
from expyriment.stimuli import Canvas
from ..pil_image.pil_image import PILImagePlotter
from .._lib.colour import Colour


class ExprimentDotArray(Canvas):

    def __init__(self, dot_array,
                 # pil_image_generator TODO better using generator
                 position=(0, 0),
                 colour_target_area=None,
                 colour_field_area=None,
                 colour_field_area_outer=None,
                 colour_center_of_mass=None,
                 colour_center_of_outer_positions=None,
                 antialiasing=True,
                 colour_background=(0, 0, 0)):
        Canvas.__init__(self, size=(0, 0), position=position)
        self.dot_array = dot_array
        self.colour_field_area = Colour(colour_field_area)
        self.colour_field_area_outer = Colour(colour_field_area_outer)
        self.colour_center_of_mass = Colour(colour_center_of_mass)
        self.colour_center_of_outer_positions = Colour(colour_center_of_outer_positions)
        self.colour_target_area = Colour(colour_target_area)
        self.colour_background = Colour(colour_background)
        self.antialiasing = antialiasing
        self._image = None

    @property
    def image(self):
        if self._image is None:
            self.make_pil_image()

        return self._image

    def make_pil_image(self):
        gen = PILImagePlotter(colour_target_area=self.colour_target_area,
                              colour_field_area=self.colour_field_area,
                              colour_field_area_outer=self.colour_field_area_outer,
                              colour_center_of_mass=self.colour_center_of_mass,
                              colour_center_of_outer_positions=self.colour_center_of_outer_positions,
                              antialiasing=self.antialiasing,
                              colour_background=self.colour_background)
        self._image = gen.plot(dot_array=self.dot_array)
        return self._image

    def _create_surface(self):
        self._size = self.image.size
        return pygame.image.frombuffer(self.image.tobytes(),
                                       self.image.size,
                                       self.image.mode)


class ExpyrimentDASequence(object):

    def __init__(self, da_sequence,
                 # pil_image_generator TODO better using generator
                 position=(0, 0),
                 colour_target_area=None,
                 antialiasing=None,
                 colour_background=(0, 0, 0),
                 make_pil_images_now=False,
                 multiprocessing=False):

        self.da_sequence = da_sequence
        self.stimuli = []
        self.position = position
        self.colour_area = colour_target_area
        self.antialiasing = antialiasing
        self.colour_background = colour_background

        for da in self.da_sequence.dot_arrays:
            stim = ExprimentDotArray(dot_array=da, position=position,
                                     colour_target_area=colour_target_area, antialiasing=antialiasing,
                                     colour_background=colour_background)
            self.stimuli.append(stim)

        if make_pil_images_now:

            if not multiprocessing:
                list(map(lambda x: x.make_pil_image(), self.stimuli))
                self._make_image_process = None
            else:
                p = Pool()

                for c, pil_im in enumerate(p.imap(ExpyrimentDASequence._make_stimuli_map_helper, self.stimuli)):
                    self.stimuli[c]._image = pil_im
                p.close()
                p.join()

    def get_stimulus_numerosity(self, number_of_dots):
        """returns stimulus with a particular numerosity"""
        try:
            return self.stimuli[self.da_sequence.numerosity_idx[number_of_dots]]
        except:
            return None

    @property
    def is_preloaded(self):
        for x in reversed(self.stimuli):
            if not x.is_preloaded:
                return False
        return True

    def preload(self, until_percent=100, time=None, do_not_return_earlier=False):
        """
        preloaded all dot_array stimuli

        Note: this will take a while!

        preload certain percent or or a time.

        """
        if until_percent > 0 and until_percent < 100:
            last = int(math.floor(until_percent * len(self.stimuli) / 100.0))
        elif until_percent == 0:
            last = 0
        else:
            last = len(self.stimuli)
        cl = Clock()

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
        return x.make_pil_image()
