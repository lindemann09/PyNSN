from __future__ import absolute_import, print_function, division

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'


import math
import pygame

from expyriment.stimuli import Canvas
from expyriment.misc import Clock
from . import pil_image


class ExprimentDotArray(Canvas):

    def __init__(self, dot_array,
                       position=(0,0),
                       colour_area=None,
                       colour_convex_hull_positions=None,
                       colour_convex_hull_dots=None,
                       colour_center_of_mass = None,
                       colour_center_of_outer_positions=None,
                       antialiasing=None,
                       colour_background=(0, 0, 0)):

        Canvas.__init__(self, size=(0, 0), position=position)
        self.dot_array = dot_array
        self.colour_convex_hull_positions  = colour_convex_hull_positions
        self.colour_convex_hull_dots  = colour_convex_hull_dots
        self.colour_center_of_mass = colour_center_of_mass
        self. colour_center_of_outer_positions = colour_center_of_outer_positions
        self.colour_area = colour_area
        self.colour_background = colour_background
        self.antialiasing = antialiasing


    def _create_surface(self):

        image = pil_image.create(dot_array=self.dot_array,
                                      colour_area=self.colour_area,
                                      colour_convex_hull_positions=self.colour_convex_hull_positions,
                                      colour_convex_hull_dots=self.colour_convex_hull_dots,
                                      colour_center_of_mass=self.colour_center_of_mass,
                                      colour_center_of_outer_positions=self.colour_center_of_outer_positions,
                                      antialiasing=self.antialiasing,
                                      colour_background=self.colour_background)

        self._size = image.size
        return pygame.image.fromstring(image.tobytes(),
                                       image.size, image.mode)


class ExpyrimentDASequence(object):

    def __init__(self, da_sequence,
                 position=(0, 0),
                 colour_area=None,
                 antialiasing=None,
                 colour_background=(0, 0, 0)):

        self.da_sequence = da_sequence
        self.stimuli = []
        for da in self.da_sequence.dot_arrrays:
            self.stimuli.append(ExprimentDotArray(dot_array=da,
                       position=position,
                       colour_area=colour_area,
                       antialiasing=antialiasing,
                       colour_background=colour_background))

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

    def preload(self, percent = 100, time=None, do_not_return_earlier=False):
        """
        preloaded all dot_array stimuli

        Note: this will take a while!

        preload certain percent or or a time.

        """

        if percent>0 and percent<100:
            last = int(math.floor(percent*len(self.stimuli)/100.0))
        elif percent==0:
            last = 0
        else:
            last = len(self.stimuli)
        cl = Clock()

        try:
            for x in self.stimuli[:last]:
                if not x.is_preloaded and (time is None or cl.time<time):
                    x.preload()
            rtn = True
        except:
            rtn = False

        if do_not_return_earlier:
            cl.wait(time - cl.time)
        return rtn


    def unload(self):
        """
        returns array of preloaded dot_array_sequence
        """

        try:
            list(map(lambda x:x.unload(), self.stimuli))
            return True
        except:
            return False
