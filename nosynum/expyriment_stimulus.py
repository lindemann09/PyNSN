from __future__ import absolute_import, print_function, division
from builtins import map, zip, filter

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


class ExpyrimentDASequence(object):

    def __init__(self, da_sequence):

        self.da_sequence = da_sequence
        self.stimuli = []

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
