
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math
from .._lib.coordinate2D import Coordinate2D


class Dot(Coordinate2D):

    def __init__(self, xy, diameter, attribute=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        xy : tuple of two numeric
        diameter : numeric
        attribute : attribute (string, optional)
        """

        Coordinate2D.__init__(self, x=xy[0], y=xy[1])
        self.diameter = diameter
        if attribute is not None and not isinstance(attribute, str):
            raise ValueError("attributes must be a string or None, not {}".format(type(attribute).__name__))
        self.attribute = attribute

    def __repr__(self):
        return "Dot(xy={}, diameter={}, attribute={})".format(self.xy,
                            self.diameter, repr(self.attribute))

    def distance(self, d):
        """Return Euclidean distance to the dot d. The function takes the
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


class Rectangle(Coordinate2D):

    def __init__(self, xy, size, attribute=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        xy : tuple
            tuple of two numeric (default=(0, 0))
        size : tuple
            tuple of two numeric (default=(0, 0))
        attribute : attribute (string)
        """

        Coordinate2D.__init__(self, x=xy[0], y=xy[1])
        if attribute is not None and not isinstance(attribute, str):
            raise ValueError("attributes must be a string or None, not {}".format(type(attribute).__name__))
        self.attribute = attribute
        self.width, self.height  = size

    def __repr__(self):
        return "Rectangle(xy={}, size={}, attribute={})".format(self.xy,
                                    self.size, repr(self.attribute))

    @property
    def left(self):
        return self.x - 0.5 * self.width

    @property
    def top(self):
        return self.y + 0.5 * self.height

    @property
    def right(self):
        return self.x + 0.5 * self.width

    @property
    def bottom(self):
        return self.y - 0.5 * self.height

    def edges(self):
        """iterator over Coordinate2D representing all four edges
        """
        yield Coordinate2D(self.left, self.top)
        yield Coordinate2D(self.right, self.top)
        yield Coordinate2D(self.right, self.bottom)
        yield Coordinate2D(self.left, self.bottom)

    @property
    def size(self):
        return self.width, self.height

    @size.setter
    def size(self, size):
        self.width, self.height  = size

    def is_point_inside_rect(self, xy):
        return (self.left <= xy[0] <= self.right and
                self.top <= xy[1] <= self.bottom)

    def overlaps_with(self, rect):
        d = self.xy_distances(rect)
        return not(d[0]>0 or d[1]>0)

    @property
    def area(self):
        return self.width * self.height

    @property
    def perimeter(self):
        return 2 * (self.width + self.height)

    @property
    def diagonal(self):
        return math.sqrt(self.width**2 + self.height**2)

    def xy_distances(self, other):
        """return distances on both axes between rectangles. 0 indicates
        overlap off edges along that dimension.
        """
        assert isinstance(other, Rectangle)
        #overlaps in x or y
        pos_dist = abs(self.x - other.x), abs(self.y - other.y)
        max_overlap_dist = (self.width + other.width)/2, (self.height + other.height)/2
        if  pos_dist[0] <= max_overlap_dist[0]:
            dx = 0
        else:
            dx = pos_dist[0] - max_overlap_dist[0]

        if  pos_dist[1] <= max_overlap_dist[1]:
            dy = 0
        else:
            dy = pos_dist[1] - max_overlap_dist[1]

        return dx, dy

    def distance(self, other):
        """Return Euclidean distance to the dot d. The function takes the
        size of the rectangles into account.

        Parameters
        ----------
        other : Rectangle

        Returns
        -------
        distance : float

        """
        dx, dy = self.xy_distances(other=other)
        return math.hypot(dx, dy)

