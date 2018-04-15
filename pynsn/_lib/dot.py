"""
Dot Array
"""
from __future__ import print_function, division
from builtins import *

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math
from .colours import convert_colour


class Coordinate2D(object):

    def __init__(self, x=0, y=0):
        self._x = x
        self._y = y
        self._pos_radius = None
        self._pos_angle = None

    @property
    def x(self):
        if self._x is None:
            self._x = self._pos_radius * math.cos(self._pos_angle)
        return self._x

    @x.setter
    def x(self, value):
        self.xy = (value, self.y)

    @property
    def y(self):
        if self._y is None:
            self._y = self._pos_radius * math.sin(self._pos_angle)
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
            self._pos_angle = math.atan2(self._y, self._x)
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
        """Return Euclidean distance to the another Coordinate. The function take the
        diameter of the points into account.

        Parameters
        ----------
        d : Dot

        Returns
        -------
        distance : float

        """

        return math.hypot(self.x - d.x, self.y - d.y)



class _Item(Coordinate2D):

    def __init__(self, x=0, y=0, colour=None, picture=None):
        """Base class for dot and rectangles

        has a colour and/or a picture
        """

        Coordinate2D.__init__(self, x=x, y=y)
        self.picture = picture
        self.colour = colour

    @property
    def colour(self):
        """RGB colour """
        return self._colour

    @colour.setter
    def colour(self, value):
        self._colour = convert_colour(value)



class Dot(_Item):

    def __init__(self, x=0, y=0, diameter=1, colour=None, picture=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        diameter : numeric (default=1)
        colour : rgb colour (default=None)
        picture : todo

        """

        _Item.__init__(self, x=x, y=y, colour=colour, picture=picture)

        self.diameter = diameter

    def __repr__(self):
        rtn = "[{0},{1}], d={2}".format(self._x, self._y, self.diameter)
        if self.colour is not None:
            rtn += ", c={0}".format(self.colour)
        return "(" + rtn + ")"


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

        return Coordinate2D.distance(self, d) - \
               ((self.diameter + d.diameter) / 2.0)

    @property
    def area(self):
        return math.pi * (self.diameter**2)/4.0

    @property
    def circumference(self):
        return math.pi * self.diameter


class Rectangle(_Item):
    def __init__(self, x=0, y=0, width=0, height=0, colour=None, picture=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        width : numeric (default=1)
        height : numeric (default=1)
        colour : rgb colour (default=None)
        picture : todo
        """

        _Item.__init__(self, x=x, y=y, colour=colour, picture=picture)

        self.height = height
        self.width = width

    @property
    def left_top(self):
        return self._x - 0.5*self.width, self._y + 0.5*self.height

    @property
    def right_bottom(self):
        return self._x + 0.5*self.width, self._y - 0.5*self.height

    @property
    def rect(self):
        rtn = list(self.left_top)
        rtn.extend(self.right_bottom)
        return rtn

    def __repr__(self):
        rtn = "[{0},{1}], w={2}, h={3}".format(self._x, self._y, self.width, self.height)
        if self.colour is not None:
            rtn += ", c={0}".format(self.colour)
        return "(" + rtn + ")"

    def distance(self, r):
        """Return Euclidean distance to the dot d. The function take the
        diameter of the points into account.

        Parameters
        ----------
        d : Dot

        Returns
        -------
        distance : float

        """

        l1, t1 = self.left_top
        r1, b1 = self.right_bottom
        l2, t2 = r.left_top
        r2, b2 = r.right_bottom

        return Coordinate2D.distance(self, d)


    def gap_xy(self, other):
        """Gap beween two rectangles seppart for the x and y axis"""

        #  overlaps in x or y:
        if abs(self.x - other.x) <= (self.width + other.width):
            dx = 0;
        else:
            dx = abs(self.x - other.x) - (self.width + other.width)
        #
        if abs(self.y - other.y) <= (self.height + other.height):
            dy = 0;
        else:
            dy = abs(self.y - other.y) - (self.h + other.height)

        return dx, dy


    @property
    def area(self):
        return self.width * self.height

    @property
    def circumference(self):
        return 2 * (self.width+self.height)

