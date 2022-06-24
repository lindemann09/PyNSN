__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math

class Coordinate2D(object):

    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __repr__(self):
        return "Coordinate2D(xy={})".format(self.xy)

    def __add__(self, other):
        return Coordinate2D(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Coordinate2D(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return Coordinate2D(self.x * other, self.y * other)

    def __div__(self, other):
        return Coordinate2D(self.x / other if other else self.x,
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
    def xy(self):
        return self.x, self.y

    @xy.setter
    def xy(self, xy_tuple):
        self.x = xy_tuple[0]
        self.y = xy_tuple[1]

    @property
    def polar_radius(self):
        return math.hypot(self.x, self.y)

    @polar_radius.setter
    def polar_radius(self, value):
        self.polar = (value, self.polar_angle)

    @property
    def polar_angle(self):
        return math.atan2(self.y, self.x)

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

        self.x = rad_ang[0] * math.cos(rad_ang[1])
        self.y = rad_ang[0] * math.sin(rad_ang[1])


    def distance(self, d):
        """Returns Euclidean distance to the another Coordinate. The function
        does not takes the size of an object into account.

        Parameters
        ----------
        d : Coordinate2D

        Returns
        -------
        distance : float

        """

        return math.hypot(self.x - d.x, self.y - d.y)
