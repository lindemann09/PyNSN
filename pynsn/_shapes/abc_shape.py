from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from typing import Any, Tuple, Union

from shapely import Polygon

from ..image._image_colours import Colour
from .._lib.geometry import Coord2DLike


class ShapeType(metaclass=ABCMeta):
    """Abstract Shape Type Class"""

    __slots__ = ("attribute", "_xy", "_polyg", "_polyg_buffered", "_buffer")

    def __init__(self, xy: Coord2DLike, attribute: Any) -> None:
        self._xy = tuple(xy)
        if len(self._xy) != 2:
            raise ValueError(
                "size has be an iterable object of two numerals (width, height)"
            )
        self.attribute = attribute
        # caches polygons
        self._polyg = None
        self._polyg_buffered = None
        self._buffer = 0  # last buffer size for poly_buffered

    @property
    def xy(self) -> Tuple:
        return self._xy

    @xy.setter
    def xy(self, val: Coord2DLike):
        self._xy = tuple(val)
        if len(self._xy) != 2:
            raise ValueError(
                "size has be an iterable object of two numerals (width, height)"
            )
        self._clear_cached_polygons()

    def polygon(self, buffer: int = 0) -> Polygon:
        """return shapely polygon of this object"""

        if buffer == 0:
            if isinstance(self._polyg, Polygon):
                return self._polyg
            else:
                self._polyg = self._make_polygon(0)
                return self._polyg
        else:
            if isinstance(self._polyg_buffered, Polygon) and self._buffer == buffer:
                return self._polyg_buffered
            else:
                self._polyg_buffered = self._make_polygon(buffer)
                self._buffer = buffer
                return self._polyg_buffered

    def _clear_cached_polygons(self):
        self._polyg = None
        self._polyg_buffered = None  # caches polygons

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

    def intersects(self, other: Union[ShapeType, Polygon]) -> bool:
        """Returns True if geometries intersect, else False

        The function wraps shapes.polygon.intersects
        """
        if isinstance(other, ShapeType):
            other = other.polygon()
        return self.polygon().intersects(other)
