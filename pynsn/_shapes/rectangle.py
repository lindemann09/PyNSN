from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Any, Iterator, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .. import _shapes
from .._lib.coordinate import Coordinate
from .abc_shape import ABCShape


class Rectangle(ABCShape):
    __slot__ = ("_size", )

    def __init__(self, xy: ArrayLike = (0, 0),
                 size: ArrayLike = (0, 0),
                 attribute: Any = None):
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
        self._size = np.empty(2)
        self.size = size  # call setter

    def __repr__(self):
        return f"Rectangle(xy={self.xy}, size={self.size}, " \
            + "attribute='{self.attribute}')"

    @property
    def width(self) -> float:
        return self._size[0]

    @property
    def height(self) -> float:
        return self._size[1]

    @property
    def area(self) -> float:
        return self.width * self.height

    @property
    def perimeter(self) -> float:
        return 2 * (self.width + self.height)

    def distance(self, other: Union[_shapes.ShapeType, Coordinate]) -> float:
        # inherited doc
        if isinstance(other, _shapes.Dot):
            dist = self.distance(_shapes.Point(xy=other.xy))
            return dist - other.diameter / 2.0

        elif isinstance(other, (_shapes.Rectangle, _shapes.Point)) \
                or other.__class__ == Coordinate:
            d_xy = self._xy_distances(other=other)
            return np.hypot(d_xy[0], d_xy[1])

        raise NotImplementedError(f"distance to {type(other)} "
                                  + "is implemented.")

    @property
    def left(self) -> float:
        return self._xy[0] - 0.5 * self._size[0]

    @property
    def top(self) -> float:
        return self._xy[1] + 0.5 * self._size[1]

    @property
    def right(self) -> float:
        return self._xy[0] + 0.5 * self._size[0]

    @property
    def bottom(self) -> float:
        return self._xy[1] - 0.5 * self._size[1]

    def get_ltrb(self) -> NDArray:
        '''returns coordinates of Left-Top and Right-Bottom edges

        Returns:
            NDArray with coordinates [[left, top],
                                      [right, botton]]
        '''
        rtn = np.array([self._xy - (self._size / 2),
                        self._xy + (self._size / 2)])
        rtn[:, 1] = np.flip(rtn[:, 1])
        return rtn

    def iter_edges(self) -> Iterator[Coordinate]:
        """iterator over all four edges

        Returns
        -------
        iterator over Coordinates or tuple (x,y)
        """

        edges = self.get_ltrb()
        # left top
        yield Coordinate(xy=edges[0, :])
        # right top
        yield Coordinate(xy=np.array((edges[0, 1], edges[1, 0])))
        # right bottom
        yield Coordinate(xy=edges[1, :])
        # left bottom
        yield Coordinate(xy=edges.diagonal())

    @property
    def size(self) -> NDArray:
        return self._size

    @size.setter
    def size(self, values: ArrayLike) -> None:
        values = np.asarray(values)
        if values.shape != (2,):
            raise ValueError(
                "size has be an iterable object with two elements (width, height)")
        self._size = values

    @property
    def proportion(self) -> float:
        """Proportion of the rectangle (width/height)"""
        return self._size[0] / self._size[1]

    @property
    def diagonal(self) -> None:
        '''size of the diagonal'''
        return np.sqrt(self._size[0] ** 2 + self._size[1] ** 2)

    def _xy_distances(self, other: Union[_shapes.ShapeType, Coordinate]) -> NDArray:
        """return distances on both axes between rectangles. 0 indicates
        overlap of edges along that dimension.
        """
        if isinstance(other, _shapes.Rectangle):
            max_overlap_dist = (self._size + other.size) / 2
        elif isinstance(other, _shapes.Point) \
                or other.__class__ == Coordinate:
            max_overlap_dist = self._size / 2
        else:
            raise NotImplementedError(f"xy_distances to {type(other)} "
                                      + "is implemented.")
        # overlaps in x or y
        pos_dist = np.abs(self._xy - other.xy) - max_overlap_dist
        pos_dist[np.where(pos_dist < 0)] = 0
        return pos_dist

    def is_inside(self, other: Union[_shapes.ShapeType, Coordinate]) -> bool:
        # inherited doc
        if isinstance(other, _shapes.Dot):
            r_other = other.diameter / 2.0
            for edge in self.iter_edges():
                # distance between other center & edge > radius
                if Coordinate.distance(other, edge) > r_other:
                    return False
            return True

        elif isinstance(other, _shapes.Rectangle):
            is_smaller = self.get_ltrb() < other.get_ltrb()
            if is_smaller[0, 0] \
                    or not is_smaller[0, 1] \
                    or not is_smaller[1, 0] \
                    or is_smaller[1, 1]:
                return False
            else:
                return True

        elif isinstance(other, _shapes.Point) or other.__class__ == Coordinate:
            return False

        raise NotImplementedError("is_inside is not "
                                  f"implemented for {type(other)}.")

# TODO doc "basic properities" is all classes (like width)
