from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from typing import Any, Union

from shapely import Polygon

from ..image._image_colours import Colour
from .._lib.geometry import Coord2DLike, Coord2D


class ShapeType(metaclass=ABCMeta):
    """Abstract Shape Type Class"""

    __slots__ = ("attribute", "_xy", "_polyg", "_polyg_buffered", "_buffer")

    def __init__(self, xy: Coord2DLike, attribute: Any) -> None:
        if len(xy) != 2:
            raise ValueError(
                "size has be an iterable object of two numerals (width, height)"
            )
        self._xy = xy
        self.attribute = attribute
        # caches polygons
        self._polyg = None
        self._polyg_buffered = None
        self._buffer = 0  # last buffer size for poly_buffered

    @property
    def xy(self) -> Coord2D:
        return self._xy  # type: ignore

    def polygon(self, buffer: int = 0) -> Polygon:
        """return shapely polygon of this object"""

        if buffer == 0:
            if not isinstance(self._polyg, Polygon):
                self._polyg = self._make_polygon(0)
            return self._polyg
        else:
            if not isinstance(self._polyg_buffered, Polygon) or self._buffer != buffer:
                self._polyg_buffered = self._make_polygon(buffer)
                self._buffer = buffer
            return self._polyg_buffered

    @property
    def colour(self) -> Colour:
        """Class instance of the attribute, if possible

        Returns
        -------
        rtn : Colour
        """
        try:
            return Colour(self.attribute)
        except TypeError:
            return Colour(None)

    @abstractmethod
    def __repr__(self) -> str:
        """"""

    @abstractmethod
    def _make_polygon(self, buffer: int = 0) -> Polygon:
        """make a polygon"""

    @property
    @abstractmethod
    def area(self) -> float:
        """area of the object/shape"""

    @property
    @abstractmethod
    def perimeter(self) -> float:
        """perimeter"""

    def intersects(self, other: Union[ShapeType, Polygon], buffer: int = 0) -> bool:
        """Returns True if geometries intersect, else False

        The function wraps shapes.polygon.intersects
        """
        if isinstance(other, ShapeType):
            other = other.polygon()
        return self.polygon(buffer).intersects(other)

    def is_inside(self, other: Union[ShapeType, Polygon], buffer: int = 0) -> bool:
        """Returns True if the shape is fully inside the other, else False

        The function wraps shapes.polygon.covered_by
        """
        if isinstance(other, ShapeType):
            other = other.polygon()
        return self.polygon(buffer).covered_by(other)

    def distance(self, other: Union[ShapeType, Polygon]) -> float:
        """Unitless distance to other geometry (float)

        The function wraps shapes.polygon.distance
        """
        if isinstance(other, ShapeType):
            other = other.polygon()
        return self.polygon().distance(other)
