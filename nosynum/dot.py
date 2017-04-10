"""
Dot Array
"""

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math

class Dot(object):

    def __init__(self, x=0, y=0, diameter=1, colour=None, picture=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        diameter : numeric (default=1)
        colour : colour (default=None)

        """

        self._x = x
        self._y = y
        self._radius = None
        self._angle = None
        self.diameter = diameter
        self.colour = colour
        self.picture = picture

    def __repr__(self):
        rtn = "[{0},{1}], d={2}".format(self._x, self._y, self.diameter)
        if self.colour is not None:
            rtn += ", c={0}".format(self.colour)
        return "(" + rtn + ")"

    @property
    def x(self):
        if self._x is None:
            self._x = self._radius * math.sin(self._angle)
        return self._x

    @x.setter
    def x(self, value):
        self.xy = (value, self.y)

    @property
    def y(self):
        if self._y is None:
            self._y = self._radius * math.cos(self._angle)
        return self._y

    @y.setter
    def y(self, value):
        self.xy = (self.x, value)

    @property
    def xy(self):
        return (self.x, self.y)

    @xy.setter
    def xy(self, xy_tuple):
        self._x = xy_tuple[0]
        self._y = xy_tuple[1]
        self._radius = None
        self._angle = None

    @property
    def radius(self):
        if self._radius is None:
            self._radius = math.hypot(self._x, self._y)
        return self._radius

    @radius.setter
    def radius(self, value):
        self.polar = (value, self.angle)

    @property
    def angle(self):
        if self._angle is None:
            self._angle = math.atan2(self._x, self._y)
        return self._angle

    @angle.setter
    def angle(self, value):
        self.polar = (self.radius, value)

    @property
    def polar(self):
        """polar coordinate (radius, angle) """
        return (self.radius, self.angle)

    @polar.setter
    def polar(self, rad_ang):
        """polar coordinate (radius, angle) """

        self._radius = rad_ang[0]
        self._angle = rad_ang[1]
        self._x = None
        self._y = None



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

        return math.hypot(self._x - d._x, self._y - d._y) - \
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
