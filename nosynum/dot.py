"""
Dot Array
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math

class Dot(object):

    def __init__(self, x=0, y=0, diameter=1, colour=None, picture=None):
        """Initialize a point

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        diameter : numeric (default=1)
        colour : colour (default=None)

        """

        self.x = x
        self.y = y
        self.diameter = diameter
        self.colour = colour
        self.picture = picture

    def __repr__(self):
        rtn = "[{0},{1}], d={2}".format(self.x, self.y, self.diameter)
        if self.colour is not None:
            rtn += ", c={0}".format(self.colour)
        return "(" + rtn + ")"

    @property
    def xy(self):
        return (self.x, self.y)

    @xy.setter
    def xy(self, xy_tuple):
        self.x = xy_tuple[0]
        self.y = xy_tuple[1]

    def distance(self, d):
        """Return Euclidean distance to the dot d. The function take the
        diameter of the points into account.

        Parameters
        ----------
        d : Dot

        Returns
        -------
        distance : float

        """

        return math.hypot(self.x - d.x, self.y - d.y) - \
              ((self.diameter + d.diameter) / 2.0)

    @property
    def radius(self):
        return float(self.diameter)/2

    @property
    def area(self):
        return math.pi * (self.diameter**2)/4

    @property
    def circumference(self):
        return math.pi * self.diameter
