from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Any, Iterator, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .. import _shapes
from .._lib.coordinate import Coordinate
from .._lib.spatial_relations import RectangleRectangle
from .abc_shape import ABCShape

# TODO doc "basic properities" is all classes (like width)


class Rectangle(ABCShape):
    __slot__ = ("_size", )

    def __init__(self, xy: ArrayLike = (0, 0),
                 size: ArrayLike = (0, 0),
                 attribute: Any = None):
        """Initialize a Rectangle

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
            dist = self.distance(Coordinate(xy=other.xy))
            return dist - other.diameter / 2.0

        elif isinstance(other, _shapes.Rectangle):
            rel = RectangleRectangle(a_xy=self._xy,
                                     a_sizes=self._size,
                                     b_xy=other.xy,
                                     b_sizes=other.size)
            return rel.distances()[0]

        elif other.__class__ == Coordinate:
            rel = RectangleRectangle(a_xy=self._xy,
                                     a_sizes=self._size,
                                     b_xy=other.xy,
                                     b_sizes=np.zeros(2))
            return rel.distances()[0]

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
        '''returns coordinates of Left-Top and Right-Bottom corners

        Returns:
            NDArray with coordinates [[left, top],
                                      [right, bottom]]
        '''
        rtn = np.array([self._xy - (self._size / 2),
                        self._xy + (self._size / 2)])
        rtn[:, 1] = np.flip(rtn[:, 1])
        return rtn

    def iter_corners(self) -> Iterator[Coordinate]:
        """iterator over all four corners

        Returns
        -------
        iterator over Coordinates or tuple (x,y)
        """

        corners = self.get_ltrb()
        # left left, top
        yield Coordinate(xy=corners[0, :])
        # right, top
        yield Coordinate(xy=np.array((corners[1, 0], corners[0, 1])))
        # right, bottom
        yield Coordinate(xy=corners[1, :])
        #  left, bottom
        yield Coordinate(xy=corners.diagonal())

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

    def is_inside(self, other: _shapes.ShapeType) -> bool:
        # FIXME remove this function
        if isinstance(other, _shapes.Dot):
            r_other = other.diameter / 2.0
            for corner in self.iter_corners():
                # distance between other center & corner > radius
                if Coordinate.distance(other, corner) > r_other:
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

    def rectangles_inside(self, xy: NDArray, sizes: NDArray) -> Union[np.bool_, NDArray[np.bool_]]:
        # inherited doc
        # FIXME not tested

        xy = np.atleast_2d(xy)
        rect_sizes2 = np.atleast_2d(sizes)/2

        # true if inside on this dimension
        compare_lb_rt = np.empty(shape=(xy.shape[0], 4))
        compare_lb_rt[:, (0, 1)] = xy - rect_sizes2 \
            > np.atleast_2d(self.xy - self.size / 2)
        compare_lb_rt[:, (2, 3)] = xy + rect_sizes2 \
            > np.atleast_2d(self.xy + self.size / 2)

        if compare_lb_rt.shape[0] > 1:
            return np.all(compare_lb_rt, axis=1)
        else:
            return np.all(compare_lb_rt)

    def dots_inside(self, xy: NDArray, diameters: NDArray) -> Union[np.bool_, NDArray[np.bool_]]:
        # inherited doc
        # FIXME not tested

        xy = np.atleast_2d(xy)
        radius = np.atleast_2d(diameters/2)
        # if true if dot inside rect on this dimension
        compare_lb_rt = np.empty(shape=(xy.shape[0], 4))
        # type: ignore
        compare_lb_rt[:, (0, 1)] = xy - \
            radius.T > np.atleast_2d(self.xy - self.size / 2)
        compare_lb_rt[:, (2, 3)] = xy + \
            radius.T > np.atleast_2d(self.xy + self.size / 2)

        if compare_lb_rt.shape[0] > 1:
            return np.all(compare_lb_rt, axis=1)
        else:
            return np.all(compare_lb_rt)
