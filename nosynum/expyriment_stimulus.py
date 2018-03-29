from __future__ import absolute_import, print_function, division
from builtins import *

__all__ = ["create", "ExpyrimentDASequence"]
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math
from expyriment.stimuli import Canvas, Circle, Line, Picture
from expyriment.misc import Clock
try:
    from .expyriment_pil_image import ExprimentPILImage
except:
    ExprimentPILImage = None
    print("Pillow (PIL) is not installed. Create_stimuli_from_pil_images will not be possible.")


def create(dot_array, area_colour=None,
           convex_hull_colour=None,
           anti_aliasing=None,
           background_stimulus_expyriment=None):

    if background_stimulus_expyriment is None:
        canvas = Canvas(size=(dot_array.definition.stimulus_area_radius * 2,) * 2)
    else:
        canvas = background_stimulus_expyriment.copy()

    if area_colour is not None:
        Circle(radius=dot_array.definition.stimulus_area_radius,
               colour=area_colour,
               anti_aliasing=anti_aliasing).plot(canvas)
    if convex_hull_colour is not None:
        # plot convey hull
        hull = dot_array.convex_hull_points
        hull = list(hull) + [hull[0]]
        last = None
        for p in hull:
            if last is not None:
                Line(start_point=last, end_point=p, line_width=2,
                     colour=convex_hull_colour,
                     anti_aliasing=anti_aliasing).plot(canvas)
            last = p

    # plot dots
    for d in dot_array.dots:
        if d.picture is not None:
            Picture(filename=d.picture,
                    position=d.xy).plot(canvas)
        else:
            Circle(radius=d.diameter//2, colour=d.colour,
                   line_width=0, position=d.xy,
                   anti_aliasing=anti_aliasing).plot(canvas)

    return canvas


class ExpyrimentDASequence(object):

    def __init__(self, da_sequence):

        self.da_sequence = da_sequence
        self.stimuli = []

    def create_stimuli_native(self, area_colour=None,
                              convex_hull_colour=None,
                              anti_aliasing=None,
                              background_stimulus_expyriment=None):
        """note: rounds dot array to intergers """

        self.stimuli = []
        for da in self.da_sequence.dot_arrays:
            da.round_dot_paramter_to_integer()
            self.stimuli.append(create(dot_array=da, area_colour=area_colour,
                                       convex_hull_colour=convex_hull_colour,
                                       anti_aliasing=anti_aliasing,
                                       background_stimulus_expyriment=background_stimulus_expyriment))

    def load_associated_images(self, position=(0,0)):
        """create stimuli from images, if images comprises filenames"""
        self.stimuli = []
        for image in self.da_sequence.images:
            self.stimuli.append(Picture(filename=image, position=position))

    def create_stimuli_from_pil_images(self, position=(0, 0), delete_pil_images_afterwards=False):
        """create stimuli from images, if images comprises pil images"""
        self.stimuli = []
        for image in self.da_sequence.images:
            self.stimuli.append(ExprimentPILImage(pil_image=image, position=position))

        if delete_pil_images_afterwards:
            self.da_sequence.images = []


    def get_stimulus_numerosity(self, number_of_dots):
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
        returns array of preloaded dot_array_sequence

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
