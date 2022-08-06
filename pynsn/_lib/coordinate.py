__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math as _math
from typing import Tuple


class Coordinate(object):

    def __init__(self, x:float =0, y:float =0):
        self.x = x
        self.y = y

    def __repr__(self):
        return "Coordinate(xy={})".format(self.xy)

    def __add__(self, other):
        return Coordinate(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Coordinate(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return Coordinate(self.x * other, self.y * other)

    def __div__(self, other):
        return Coordinate(self.x / other if other else self.x,
                          self.y / other if other else self.y)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        return self

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        return self

    def __imul__(self, other):
        self.x *= other
        self.y *= other
        return self

    def __idiv__(self, other):
        self.x /= other
        self.y /= other
        return self

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return self.x != other.x or self.y != other.y

    @property
    def xy(self) -> Tuple[float, float]:
        return self.x, self.y

    @xy.setter
    def xy(self, xy_tuple):
        self.x = xy_tuple[0]
        self.y = xy_tuple[1]

    @property
    def polar_radius(self):
        return _math.hypot(self.x, self.y)

    @polar_radius.setter
    def polar_radius(self, value):
        self.polar = (value, self.polar_angle)

    @property
    def polar_angle(self):
        return _math.atan2(self.y, self.x)

    @polar_angle.setter
    def polar_angle(self, value):
        self.polar = (self.polar_radius, value)

    @property
    def polar(self):
        """polar coordinate (radius, pos_angle) """
        return self.polar_radius, self.polar_angle

    @polar.setter
    def polar(self, rad_ang):
        """polar coordinate (radius, angle) """

        self.x = rad_ang[0] * _math.cos(rad_ang[1])
        self.y = rad_ang[0] * _math.sin(rad_ang[1])

    def distance(self, other):
        """Euclidean distance to the another Coordinate."""
        return _math.hypot(self.x - other.x, self.y - other.y)

# TODO typing
