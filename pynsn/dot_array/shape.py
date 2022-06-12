
__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math
from ..lib.geometry import lines_intersect, distance_between_edge_and_point
from ..lib.coordinate2D import Coordinate2D


class Dot(Coordinate2D):

    def __init__(self, x=0, y=0, diameter=1, attribute=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        diameter : numeric (default=1)
        attribute : attribute (string)
        """

        Coordinate2D.__init__(self, x=x, y=y)
        self.diameter = diameter
        if attribute is not None and not isinstance(attribute, str):
            raise TypeError("attributes must be a string or None, not {}".format(type(attribute).__name__))
        self.attribute = attribute

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

    def __init__(self, center_x=0, center_y=0, width=0, height=0, attribute=None):
        """Initialize a point

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Parameters
        ----------
        x : numeric (default=0)
        y : numeric (default=0)
        width : numeric (default=1)
        height : numeric (default=1)
        attribute : attribute (string)
        """

        Coordinate2D.__init__(self, x=center_x, y=center_y)
        if attribute is not None and not isinstance(attribute, str):
            raise TypeError("attributes must be a string or None, not {}".format(type(attribute).__name__))
        self.attribute = attribute
        self.height = height
        self.width = width

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
        """
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

