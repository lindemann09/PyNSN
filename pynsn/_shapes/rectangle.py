from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from typing import Any, List, Tuple, Optional
from math import sqrt

from shapely import Polygon

from .abc_shape import ShapeType
from .._lib.geometry import Coord2DLike, Coord2D


class Rectangle(ShapeType):
    __slots__ = ("_geo", "_size")

    QUAD_SEGS = 16  #  line segments used to approximate buffer

    def __init__(
        self,
        xy: Coord2DLike = (0, 0),
        size: Coord2DLike = (0, 0),
        attribute: Any = None,
    ):
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
        self._size = tuple(size)
        if len(self._size) != 2:
            raise ValueError("size has be an tuple of two numerals (width, height)")

    def variant(
        self,
        xy: Optional[Coord2DLike] = None,
        size: Optional[Coord2DLike] = None,
        attribute: Optional[Any] = None,
    ) -> Rectangle:
        """return a variant of this Rectangle"""
        if xy is None:
            xy = self.xy
        if size is None:
            size = self.size
        if attribute is None:
            attribute = self.attribute
        return Rectangle(xy, size, attribute)

    def _make_polygon(self, buffer: int = 0) -> Polygon:
        l = self._xy[0] - self._size[0] / 2
        r = self._xy[0] + self._size[0] / 2
        t = self._xy[1] + self._size[1] / 2
        b = self._xy[1] - self._size[1] / 2
        rtn = Polygon(((l, b), (l, t), (r, t), (r, b)))
        if buffer == 0:
            return rtn
        else:
            return rtn.buffer(buffer, quad_segs=Rectangle.QUAD_SEGS)

    def __repr__(self):
        return (
            f"Rectangle(xy={self._xy}, size={self._size}, "
            + f"attribute='{self.attribute}')"
        )

    @property
    def size(self) -> Coord2D:
        return self._size  # type: ignore

    @property
    def width(self) -> float:
        return self._size[0]

    @property
    def height(self) -> float:
        return self._size[1]

    @property
    def area(self) -> float:
        return self._size[0] * self._size[1]

    @property
    def perimeter(self) -> float:
        return 2 * (self._size[0] + self._size[1])

    @property
    def left_top(self) -> Tuple:
        return self.polygon().exterior.coords[1]

    @property
    def right_bottom(self) -> Tuple:
        return self.polygon().exterior.coords[3]

    @property
    def edges(self) -> List:
        """returns coordinates of all four edges clockwise starting with Left-Top

        Returns:
            2DArray with xy-coordinates
        """
        return self.polygon().exterior.coords[1:]

    @property
    def proportion(self) -> float:
        """Proportion of the rectangle (width/height)"""
        return self._size[0] / self._size[1]

    @property
    def diagonal(self) -> float:
        """size of the diagonal"""
        return sqrt(self._size[0] ** 2 + self._size[1] ** 2)
