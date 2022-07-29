from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import math
from typing import Iterator
from .abc_shape import ABCShape
from .._lib.coordinate import Coordinate
from .. import _shapes


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

    def distance(self, other: _shapes.ShapeType) -> float:
        # inherited doc
        if isinstance(other, _shapes.Dot):
            dist = self.distance(_shapes.Point(xy=other.xy))
            return dist - other.diameter / 2.0

        elif isinstance(other, (_shapes.Rectangle, _shapes.Point)):
            d_xy = self.xy_distances(other=other)
            return math.hypot(d_xy[0], d_xy[1])  # TODO distance to point

        raise NotImplementedError(f"distance to {type(other)} "
                                  + "is implemented.")

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
        if isinstance(other, _shapes.Rectangle):
            max_overlap_dist = (self.width + other.width) / 2,\
                (self.height + other.height) / 2
        elif isinstance(other, _shapes.Point):
            max_overlap_dist = self.width / 2, self.height / 2
        else:
            raise NotImplementedError(f"xy_distances to {type(other)} "
                                      + "is implemented.")

        # overlaps in x or y
        pos_dist = abs(self.x - other.x), abs(self.y - other.y)
        if pos_dist[0] <= max_overlap_dist[0]:
            dx = 0
        else:
            dx = pos_dist[0] - max_overlap_dist[0]

        if pos_dist[1] <= max_overlap_dist[1]:
            dy = 0
        else:
            dy = pos_dist[1] - max_overlap_dist[1]

        return dx, dy

    def is_inside(self, other: _shapes.ShapeType) -> bool:
        # inherited doc
        if isinstance(other, _shapes.Dot):
            pass

        elif isinstance(other, _shapes.Rectangle):
            pass

        elif isinstance(other, _shapes.Point):
            return False

        raise NotImplementedError("is_inside is not "
                                  "implemented for {}.".format(type(other)))

#    def is_point_inside_rect(self, xy):
#        return (other.left <= self.xy[0] <= other.right and
#                    other.top <= self.xy[1] <= other.bottom)
