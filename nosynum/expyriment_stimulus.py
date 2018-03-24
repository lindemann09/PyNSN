from __future__ import absolute_import, print_function, division
from builtins import *

__all__ = ["create", "ExpyrimentDASequence"]
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

from expyriment.stimuli import Canvas, Circle, Line, Picture

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


class ExpyrimentDASequence():

    def __init__(self, da_sequence):

        self.da_sequence = da_sequence
        self.stimuli = []

    def create_stimuli(self, area_colour=None,
                       convex_hull_colour=None,
                       anti_aliasing=None,
                       background_stimulus_expyriment=None):
        self.stimuli = []
        for da in self.da_sequence.dot_arrays:
            self.stimuli.append(create(dot_array=da, area_colour=area_colour,
                                       convex_hull_colour=convex_hull_colour,
                                       anti_aliasing=anti_aliasing,
                                       background_stimulus_expyriment=background_stimulus_expyriment))

    def load_associated_images(self, position=(0,0)):
        """e.g. if created with pil images"""
        self.stimuli = []
        for image in self.da_sequence.images:
            self.stimuli.append(Picture(filename=image, position=position))

    def get_stimulus_numerosity(self, number_of_dots):
        try:
            return self.stimuli[self.da_sequence.numerosity_idx[number_of_dots]]
        except:
            return None

    def preload(self):
        """
        returns array of preloaded dot_array_sequence
        """

        try:
            list(map(lambda x: x.preload(), self.stimuli))
            return True
        except:
            return False

    def unload(self):
        """
        returns array of preloaded dot_array_sequence
        """

        try:
            list(map(lambda x:x.unload(), self.stimuli))
            return True
        except:
            return False
