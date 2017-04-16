"""
Dot Array
"""
from __future__ import absolute_import, print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math

## TODO speed up by implementing a dot list class in numpy

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
        self._pos_radius = None
        self._pos_angle = None
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
            self._x = self._pos_radius * math.sin(self._pos_angle)
        return self._x

    @x.setter
    def x(self, value):
        self.xy = (value, self.y)

    @property
    def y(self):
        if self._y is None:
            self._y = self._pos_radius * math.cos(self._pos_angle)
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
        self._pos_radius = None
        self._pos_angle = None

    @property
    def pos_radius(self):
        if self._pos_radius is None:
            self._pos_radius = math.hypot(self._x, self._y)
        return self._pos_radius

    @pos_radius.setter
    def pos_radius(self, value):
        self.polar = (value, self.pos_angle)

    @property
    def pos_angle(self):
        if self._pos_angle is None:
            self._pos_angle = math.atan2(self._x, self._y)
        return self._pos_angle

    @pos_angle.setter
    def pos_angle(self, value):
        self.polar = (self.pos_radius, value)

    @property
    def polar(self):
        """polar coordinate (radius, pos_angle) """
        return (self.pos_radius, self.pos_angle)

    @polar.setter
    def polar(self, rad_ang):
        """polar coordinate (radius, angle) """

        self._pos_radius = rad_ang[0]
        self._pos_angle = rad_ang[1]
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

        return math.hypot(self.x - d.x, self.y - d.y) - \
               ((self.diameter + d.diameter) / 2.0)

    @property
    def area(self):
        return math.pi * (self.diameter**2)/4

    @property
    def circumference(self):
        return math.pi * self.diameter

