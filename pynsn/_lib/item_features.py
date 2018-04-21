"""
Dot Array
"""
from __future__ import print_function, division, unicode_literals
from builtins import map, zip, range, filter

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import numpy as np
from .colour import Colour


class ItemFeatures(object):

    def __init__(self, colour=None, picture=None):
        self.colour = Colour(colour)
        self.picture = picture

    def __repr__(self):
        return 'ItemFeatures(colour={}, picture={})'.format(
                    self.colour, self.picture)


class ItemFeaturesList(object):

    def __init__(self, colours=None, pictures=None):
        """If all ItemFeatures should None (but not empty) init with ItemFeatues(colours=[None]),
        otherwise the ItemFeatrues will be empty"""

        self.colours = np.array([])
        self.pictures = np.array([])
        if colours is  not None or pictures is not None:
            self.append(colours=colours, pictures=pictures)

    def __iter__(self):
        return map(lambda x: ItemFeatures(colour=x[0], picture=x[1]),
                                zip(self.colours, self.pictures))

    def __getitem__(self, key):
        return ItemFeatures(colour=self.colours[key],
                            picture=self.pictures[key])

    @property
    def length(self):
        return len(self.colours)

    def clear(self):
        self.colours = np.array([])
        self.pictures = np.array([])

    def delete(self, index):
        self.colours = np.delete(self.colours, index)
        self.pictures = np.delete(self.pictures, index)

    def repeat(self, repeats):
        """:return repeated features, equivalent to numpy.repeat"""

        rtn = ItemFeaturesList()
        rtn.colours = self.colours.repeat(repeats)
        rtn.pictures = self.pictures.repeat(repeats)
        return rtn

    def copy(self):
        rtn = ItemFeaturesList()
        rtn.colours = self.colours.copy()
        rtn.pictures = self.pictures.copy()
        return rtn

    def append_features(self, features):
        """features: None, ItemFeature, ItemFeatureList
        """

        if isinstance(features, ItemFeatures):
            self.append(colours=features.colour,
                        pictures=features.picture)
        else:
            self.append(colours=features.colours,
                    pictures=features.pictures)

    def append(self, colours, pictures):

        # ensure good numpy array or None
        pictures = numpy_vector(pictures)
        colours = numpy_vector(colours)

        if (colours is not None and pictures is not None and len(colours) != len(pictures)):
            raise RuntimeError(u"Bad shaped data: " + u"colour and picture have not the same length.")

        if colours is None and pictures is None:
            colours = [None]
            pictures = [None]
        elif pictures is None:
            pictures = [None] * len(colours)
        elif colours is None:
            colours = [None] * len(pictures)

        self.pictures = np.append(self.pictures, pictures)  # TODO check picture exist
        for c in colours:
            self.colours = np.append(self.colours, Colour(c))

    def change(self, colour=None, picture=None, indices=None):
        """allows usung color names and rgb array, since it
        converts colour """

        if isinstance(indices, int):
            indices = [indices]
        elif indices is None:
            indices = range(len(self.colours))

        if colour is not None:
            self.colours[indices] = Colour(colour)
        if picture is not None:
            self.pictures[indices] = picture



################### helper functions ###########################

def numpy_vector(x):
    """helper function:
    make an numpy vector from any element (list, arrays, and single data (str, numeric))
    nut None will not be procesed and returns None"""

    if x is None:
        return None

    x = np.array(x)
    if x.ndim == 1:
        return x
    elif x.ndim == 0:
        return x.reshape(1)  # if one element only, make a array with one element
    else:
        return x.flatten()