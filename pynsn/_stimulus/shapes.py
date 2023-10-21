from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from os import path
from typing import Any, Tuple, Union
from copy import deepcopy

import shapely
from shapely import Point, Polygon

from .._lib.colour import Colour
from .._lib.geometry import Coord2DLike


class AbstractShape(metaclass=ABCMeta):
    """Abstract Shape Type Class"""

    def __init__(self, xy: Coord2DLike, attribute: Any) -> None:
        if len(xy) != 2:
            raise ValueError("xy has be an list of two numerals (x, y)")
        self._xy = tuple(xy)
        self._attribute = attribute
        self._polygon = None

    @property
    def xy(self) -> Tuple[float, float]:
        return self._xy  # type: ignore

    @xy.setter
    def xy(self, val: Coord2DLike):
        if len(val) != 2:
            raise ValueError("xy has be an list of two numerals (x, y)")

        xy = tuple(val)
        if self._polygon is not None:
            move_xy = (xy[0] - self._xy[0], xy[1] - self._xy[1])
            self._polygon = shapely.transform(self._polygon,
                                              lambda x: x + move_xy)
        self._xy = xy

    @property
    def attribute(self) -> Any:
        return self._attribute

    @attribute.setter
    def attribute(self, val: Any):
        self._attribute = val

    @property
    def polygon(self) -> Polygon:
        if not isinstance(self._polygon, Polygon):
            self._polygon = self.make_polygon(0)
            shapely.prepare(self._polygon)  # FIXME required?
        return self._polygon

    def delete_polygon(self):
        self._polygon = None

    def set_polygon(self, polygon: Polygon) -> None:
        """This method should usually not been used, unless you know what you do.

        Polygon is automatically created, if required. Please use the property
        `polygon.
        """
        self._polygon = polygon

    @property
    def colour(self) -> Colour:
        """Class instance of the attribute, if possible

        Returns
        -------
        rtn : Colour
        """
        if isinstance(self._attribute, Colour):
            return self._attribute
        try:
            return Colour(self._attribute)
        except TypeError:
            return Colour(None)

    @abstractmethod
    def __repr__(self) -> str:
        """"""

    @abstractmethod
    def make_polygon(self, buffer) -> Polygon:
        """make a polygon"""

    @abstractmethod
    def copy(self) -> Polygon:
        """returns a copy of the object

        cached polygon will not be copied and recreated with the next call
        """

    def move(self, move_xy: Coord2DLike):
        """moves the xy position of the shape
        """

        if len(move_xy) != 2:
            raise ValueError("move_xy has be an list of two numerals (x, y)")
        self._xy = (self._xy[0] + move_xy[0],
                    self._xy[1] + move_xy[1])
        if self._polygon is not None:
            self._polygon = shapely.transform(self._polygon,
                                              lambda x: x + move_xy)


class Dot(AbstractShape):
    QUAD_SEGS = 32  # line segments used to approximate dot

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

    @property
    def diameter(self) -> float:
        return self._diameter

    @diameter.setter
    def diameter(self, val: float):
        self._diameter = val
        self._polygon = None

    def copy(self) -> Dot:
        # inherited doc
        return deepcopy(self)

    def make_polygon(self, buffer: int = 0) -> Polygon:
        # efficient buffer (and not buffer of dot polygon)
        bf = self._diameter / 2 + buffer
        return Point(self._xy).buffer(bf, quad_segs=Dot.QUAD_SEGS)

    def __repr__(self):
        return (
            f"Dot(xy={self._xy}, diameter={self._diameter}, "
            + f"attribute = {self._attribute})"
        )


class Rectangle(AbstractShape):
    QUAD_SEGS = 16  # line segments used to approximate buffer

    def __init__(self,
                 xy: Coord2DLike = (0, 0),
                 size: Coord2DLike = (0, 0),
                 attribute: Any = None,
                 ):
        super().__init__(xy, attribute)
        if len(size) != 2:
            raise ValueError(
                "size has be an list of two numerals (width, height)")
        self._size = tuple(size)

    @property
    def size(self) -> Tuple[float, float]:
        return self._size  # type: ignore

    @size.setter
    def size(self, val: Coord2DLike):
        if len(val) != 2:
            raise ValueError(
                "size has be an list of two numerals (width, height)")
        self._size = tuple(val)
        self._polygon = None

    def make_polygon(self, buffer: int = 0) -> Polygon:
        l = self._xy[0] - self._size[0] / 2
        r = self._xy[0] + self._size[0] / 2
        t = self._xy[1] + self._size[1] / 2
        b = self._xy[1] - self._size[1] / 2
        poly = Polygon(((l, b), (l, t), (r, t), (r, b)))
        if buffer != 0:
            poly = poly.buffer(buffer, quad_segs=Rectangle.QUAD_SEGS)
        return poly

    @property
    def width(self) -> float:
        return self._size[0]

    @property
    def height(self) -> float:
        return self._size[1]

    @property
    def proportion(self) -> float:
        """Proportion of the rectangle (width/height)"""
        return self._size[0] / self._size[1]

    @property
    def left(self) -> float:
        return self._xy[0] - 0.5 * self._size[0]

    @property
    def right(self) -> float:
        return self._xy[0] + 0.5 * self._size[0]

    @property
    def top(self) -> float:
        return self._xy[1] + 0.5 * self._size[1]

    @property
    def bottom(self) -> float:
        return self._xy[1] - 0.5 * self._size[1]

    def __repr__(self):
        return (f"Rectangle(xy={self._xy}, size={self._size}, "
                + f"attribute='{self._attribute}')")

    def copy(self) -> Rectangle:
        return deepcopy(self)


class Picture(Rectangle):
    ATTR_PREFIX = "p:"

    def __init__(
        self, xy: Coord2DLike = (0, 0), size: Coord2DLike = (0, 0), filename: str = ""
    ) -> None:
        """Initialize a Picture

        Rectangle can also consist of a picture

        Parameters
        ----------
        xy : tuple
            tuple of two numeric (default=(0, 0))
        size : tuple
            tuple of two numeric (default=(0, 0))
        filename : str
            path the picture file
        """

        super().__init__(xy=xy, size=size, attribute=Picture.ATTR_PREFIX + filename)

    def __repr__(self):
        return (f"Picture(xy={self._xy}, size={self._size}, " +
                f"filename='{self.filename}')")

    @property
    def filename(self) -> str:
        return self._attribute[len(Picture.ATTR_PREFIX):]

    def file_exists(self) -> bool:
        """Checks if the file exists"""
        return path.isfile(self.filename)

    @staticmethod
    def extract_filename(txt: str) -> Union[str, None]:
        """Check if text string is a picture attribute and returns the filename.
        Otherwise, returns None

        Args:
            txt: string to be checked
        """
        if isinstance(txt, str) and txt.startswith(Picture.ATTR_PREFIX):
            return txt[len(Picture.ATTR_PREFIX):]
        else:
            return None

    def copy(self) -> Picture:
        return deepcopy(self)
