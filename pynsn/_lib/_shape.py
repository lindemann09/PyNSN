
__author__ = 'Oliver Lindemann <oliver.lindemann@cognitive-psychology.eu>'

import math
from ._geometry import lines_intersect, distance_between_edge_and_point
from ._item_attributes import ItemAttributes

class _Coordinate2D(object):

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
        d : _Coordinate2D

        Returns
        -------
        distance : float

        """

        return math.hypot(self.x - d.x, self.y - d.y)


class Dot(_Coordinate2D):  # TODO becomes maybe an item

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

        _Coordinate2D.__init__(self, x=x, y=y)
        self.diameter = diameter
        if attributes is None:
            self.attributes = ItemAttributes()
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

        return _Coordinate2D.distance(self, d) - \
               ((self.diameter + d.diameter) / 2.0)

    @property
    def area(self):
        return math.pi * (self.diameter ** 2) / 4.0

    @property
    def perimeter(self):
        return math.pi * self.diameter


class Rectangle(_Coordinate2D):

    def __init__(self, center_x=0, center_y=0, width=0, height=0, attributes=None):
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

        _Coordinate2D.__init__(self, x=center_x, y=center_y)
        if attributes is None:
            self.attributes = ItemAttributes()
        else:
            self.attributes = attributes

        self.height = height
        self.width = width

    @property
    def left(self):
        return self._x - 0.5 * self.width

    @property
    def top(self):
        return self._y + 0.5 * self.height

    @property
    def right(self):
        return self._x + 0.5 * self.width

    @property
    def bottom(self):
        return self._y - 0.5 * self.height


    def iter_edges(self):
        yield self.left, self.top
        yield self.right, self.top
        yield self.right, self.bottom
        yield self.left, self.bottom

    def is_point_inside_rect(self, xy):
        return (self.left <= xy[0] <= self.right and
                self.top <= xy[1] <= self.bottom)

    def overlaps_with(self, rect):
        for corner in rect.iter_edges():
            if self.is_point_inside_rect(corner):
                return True
        for corner in self.iter_edges():
            if rect.is_point_inside_rect(corner):
                return True
        return False

    def distance(self, rect):
        """Return Euclidean distance to other rect. The function take the
        diameter of the points into account.

        Parameters
        ----------
        d : Dot

        Returns
        -------
        distance : float or -1 if overlapping

        """

        # 1. see if they overlap
        if self.overlaps_with(rect):
            return 0

        # 2. draw a line between rectangles
        line = (self.xy, rect.xy)

        # 3. find the two edges that intersect the line
        edge1 = None
        edge2 = None
        for edge in self.iter_edges():
            if lines_intersect(edge, line):
                edge1 = edge
                break
        for edge in rect.iter_edges():
            if lines_intersect(edge, line):
                edge2 = edge
                break
        assert edge1
        assert edge2

        # 4. find shortest distance between these two edges
        distances = [
            distance_between_edge_and_point(edge1, edge2[0]),
            distance_between_edge_and_point(edge1, edge2[1]),
            distance_between_edge_and_point(edge2, edge1[0]),
            distance_between_edge_and_point(edge2, edge1[1]),
        ]

        return min(distances)

    @property
    def area(self):
        return self.width * self.height

    @property
    def perimeter(self):
        return 2 * (self.width + self.height)

