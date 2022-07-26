__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math
from typing import Iterator
from .abc_shape import ABCShape
from .._lib.coordinate import Coordinate


class Rectangle(ABCShape):

    def __init__(self, xy, size, attribute=None):
        """Initialize a Rectangle

        Handles polar and cartesian representation (optimised processing, i.e.,
        conversions between coordinates systems will be done only once if needed)

        Rectangle can also consist of a picture

        Parameters
        ----------
        xy : tuple
            tuple of two numeric (default=(0, 0))
        size : tuple
            tuple of two numeric (default=(0, 0))
        attribute : attribute (string) or PictureFile
        """

        super().__init__(xy, attribute)
        self.width, self.height = size

    def __repr__(self):
        return "Rectangle(xy={}, size={}, attribute='{}')".format(self.xy,
                                    self.size, self.attribute)

    @property
    def area(self):
        return self.width * self.height

    @property
    def perimeter(self):
        return 2 * (self.width + self.height)

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

    def iter_edges(self) -> Iterator[Coordinate]:
        """iterator over all four edges

        Returns
        -------
        iterator over Coordinates or tuple (x,y)
        """

        yield Coordinate(self.left, self.top)
        yield Coordinate(self.right, self.top)
        yield Coordinate(self.right, self.bottom)
        yield Coordinate(self.left, self.bottom)

    @property
    def size(self):
        return self.width, self.height

    @size.setter
    def size(self, size):
        self.width, self.height = size

    def is_point_inside_rect(self, xy):
        return (self.left <= xy[0] <= self.right and
                self.top <= xy[1] <= self.bottom)

    def overlaps_with(self, rect):
        d = self.xy_distances(rect)
        return not (d[0] > 0 or d[1] > 0)

    @property
    def proportion(self):
        """Proportion of the rectangle (width/height)"""
        return self.width / self.height

    @property
    def diagonal(self):
        return math.sqrt(self.width ** 2 + self.height ** 2)

    def xy_distances(self, other):
        """return distances on both axes between rectangles. 0 indicates
        overlap off edges along that dimension.
        """
        assert isinstance(other, Rectangle)
        # overlaps in x or y
        pos_dist = abs(self.x - other.x), abs(self.y - other.y)
        max_overlap_dist = (self.width + other.width) / 2, (self.height + other.height) / 2
        if pos_dist[0] <= max_overlap_dist[0]:
            dx = 0
        else:
            dx = pos_dist[0] - max_overlap_dist[0]

        if pos_dist[1] <= max_overlap_dist[1]:
            dy = 0
        else:
            dy = pos_dist[1] - max_overlap_dist[1]

        return dx, dy
