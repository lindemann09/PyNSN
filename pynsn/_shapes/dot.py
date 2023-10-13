from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from typing import Any, Optional

from math import pi
from shapely import Polygon, Point

from .abc_shape import ShapeType
from .._lib.geometry import Coord2DLike


class Dot(ShapeType):
    __slots__ = "_diameter"

    QUAD_SEGS = 32  #  line segments used to approximate dot

    def __init__(
        self, xy: Coord2DLike = (0, 0), diameter: float = 0, attribute: Any = None
    ) -> None:
        """Initialize a dot

        Parameters
        ----------
        xy : tuple of two numeric (default=(0,0))
        diameter : numeric (default=(0,0))
        attribute : attribute (string, optional)
        """
        super().__init__(xy, attribute)
        self._diameter = diameter

    def variant(
        self,
        xy: Optional[Coord2DLike] = None,
        diameter: Optional[float] = None,
        attribute: Optional[Any] = None,
    ) -> Dot:
        """return a variant of this Rectangle"""
        if xy is None:
            xy = self.xy
        if diameter is None:
            diameter = self.diameter
        if attribute is None:
            attribute = self.attribute
        return Dot(xy, diameter, attribute)

    def _make_polygon(self, buffer: int = 0) -> Polygon:
        bf = self.diameter / 2 + buffer
        return Point(self._xy).buffer(bf, quad_segs=Dot.QUAD_SEGS)

    def __repr__(self):
        return (
            f"Dot(xy={self.xy}, diameter={self._diameter}, "
            + f"attribute = {self.attribute})"
        )

    @property
    def diameter(self) -> float:
        return self._diameter

    @property
    def area(self) -> float:
        return pi * (self._diameter**2) / 4.0

    @property
    def perimeter(self) -> float:
        return pi * self._diameter
