"""
Dot Array
"""
from __future__ import print_function, division
from builtins import map, range, zip, filter

__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math
from .item_attributes import ItemAttributes, ItemAttributesList


def is_all_larger(vector, standard=0):
    return sum(map(lambda x: x > standard, vector))==len(vector)

def is_all_smaller(vector, standard=0):
    return sum(map(lambda x: x < standard, vector))==len(vector)

class Coordinate2D(object):

    def __init__(self, x=0, y=0):
        self._x = x
        self._y = y
        self._polar_radius = None
        self._polar_angle = None

    @property
    def x(self):
        if self._x is None:
            self._x = self._polar_radius * math.cos(self._polar_angle)
        return self._x

    @x.setter
    def x(self, value):
        self.xy = (value, self.y)

    @property
    def y(self):
        if self._y is None:
            self._y = self._polar_radius * math.sin(self._polar_angle)
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
        self._polar_radius = None
        self._polar_angle = None

    @property
    def polar_radius(self):
        if self._polar_radius is None:
            self._polar_radius = math.hypot(self._x, self._y)
        return self._polar_radius

    @polar_radius.setter
    def polar_radius(self, value):
        self.polar = (value, self.polar_angle)

    @property
    def polar_angle(self):
        if self._polar_angle is None:
            self._polar_angle = math.atan2(self._y, self._x)
        return self._polar_angle

    @polar_angle.setter
    def polar_angle(self, value):
        self.polar = (self.polar_radius, value)

    @property
    def polar(self):
        """polar coordinate (radius, pos_angle) """
        return (self.polar_radius, self.polar_angle)

    @polar.setter
    def polar(self, rad_ang):
        """polar coordinate (radius, angle) """

        self._polar_radius = rad_ang[0]
        self._polar_angle = rad_ang[1]
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


class Dot(Coordinate2D):  # TODO becomes maybe an item

    def __init__(self, x=0, y=0, diameter=1, attributes=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        diameter : numeric (default=1)
        attributes : ItemAttributes
        """

        Coordinate2D.__init__(self, x=x, y=y)
        self.diameter = diameter
        if attributes is None:
            self.attributes = ItemAttributes(colour=None, picture=None)
        elif not isinstance(attributes, ItemAttributes):
            raise TypeError("features must be a ItemFeatures, not {}".format(type(attributes).__name__))
        else:
            self.attributes = attributes

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
        return math.pi * (self.diameter ** 2) / 4.0

    @property
    def perimeter(self):
        return math.pi * self.diameter


class Rectangle(Coordinate2D):  # todo

    def __init__(self, x=0, y=0, width=0, height=0, features=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        width : numeric (default=1)
        height : numeric (default=1)
        """

        Coordinate2D.__init__(self, x=x, y=y)
        if features is None:
            self.features = ItemAttributesList(colour=None)
        else:
            self.features = features

        self.height = height
        self.width = width

    @property
    def left_top(self):
        return self._x - 0.5 * self.width, self._y + 0.5 * self.height

    @property
    def right_bottom(self):
        return self._x + 0.5 * self.width, self._y - 0.5 * self.height

    @property
    def rect(self):
        rtn = list(self.left_top)
        rtn.extend(self.right_bottom)
        return rtn

    def distance(self, d):
        """Return Euclidean distance to the dot d. The function take the
        diameter of the points into account.

        Parameters
        ----------
        d : Dot

        Returns
        -------
        distance : float or -1 if overlapping

        """

        lt1 = self.left_top
        rb1 = self.right_bottom
        lt2 = d.left_top
        rb2 = d.right_bottom

        x_diff = [lt1[0]-lt2[0], lt1[0]-r2, r1-l2, r1-r2]
        y_diff = [t1-t2, t1-b2, b1-t2, b1-b2]
        if is_all_smaller(x_diff):
            is_right = True # other is to the right
        elif is_all_larger(x_diff):
            is_right = False
        else:
            is_right = None # is at x overlapping
        if is_all_smaller(y_diff):
            is_bottom = True # other is to the right
        elif is_all_larger(x_diff):
            is_bottom = False
        else:
            is_bottom = None # is overlapping at y

        if is_bottom is None and is_right is None:
            return -1 # overlapping or touching rects

        if is_right and is_bottom:
            a


        return Coordinate2D.distance(self, d)  # TODO check me

    def gap_xy(self, other):
        """Gap beween two rectangles seppart for the x and y axis"""


        #  overlaps in x or y:
        if abs(self.x - other.x) <= (self.width + other.width):
            dx = 0
        else:
            dx = abs(self.x - other.x) - (self.width + other.width)
        #
        if abs(self.y - other.y) <= (self.height + other.height):
            dy = 0
        else:
            dy = abs(self.y - other.y) - (self.h + other.height)

        return dx, dy

    @property
    def area(self):
        return self.width * self.height

    @property
    def perimeter(self):
        return 2 * (self.width + self.height)
