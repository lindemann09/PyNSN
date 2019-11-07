"""
Dot Array
"""
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import numpy as np

from ._colour import Colour

class ItemAttributes(object):

    def __init__(self, colour=None, picture=None):
        self.colour = colour
        self.picture = picture

    @property
    def colour(self):
        return self._colour

    @colour.setter
    def colour(self, value):
        self._colour = Colour(value)

    def __str__(self):
        return str(self.as_dict())

    def as_dict(self):
        return {"colour" : str(self.colour),
                "picture" : self.picture}

    def read_from_dict(self, dict):
        try:
            self.colour = dict["colour"]
        except:
            self.colour = None

        try:
            self.picture = dict["picture"]
        except:
            self.picture = None

    def is_different(self, other):
        """Returns True if at least on attribute, that is defined in both,
        is different. Thus, None defined attributes will be ignored"""

        assert isinstance(other, ItemAttributes)

        if self.colour != other.colour and self.colour is not None and \
                other.colour is not None:
            return True

        if self.picture != other.picture and self.picture is not None and \
                other.picture is not None:
            return True

        return False

    @staticmethod
    def check_type(data):
        """Checks if data a ItemAttributes ora list of ItemAttributes and
         raises a ValueError if not"""

        if isinstance(data, ItemAttributes):
            return True
        else:
            invalid_type_error = "Invalid type: {}.\n attributes has to be a " + \
                                 "ItemAttribute of a list of ItemAttributes."
            if isinstance(data, (list, tuple)):
                for d in data:
                    if not isinstance(d, ItemAttributes):
                        raise ValueError(invalid_type_error.format(type(d)))
                return True

            raise ValueError(invalid_type_error.format(type(data)))
