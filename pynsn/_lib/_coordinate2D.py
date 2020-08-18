__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math

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
        self._x = value
        self._polar_radius = None
        self._polar_angle = None

    @property
    def y(self):
        if self._y is None:
            self._y = self._polar_radius * math.sin(self._polar_angle)
        return self._y

    @y.setter
    def y(self, value):
        self._y = value
        self._polar_radius = None
        self._polar_angle = None


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
        d : Coordinate2D

        Returns
        -------
        distance : float

        """

        return math.hypot(self.x - d.x, self.y - d.y)


