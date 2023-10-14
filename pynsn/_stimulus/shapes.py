from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from os import path
from typing import Any, Optional, Tuple, Union

from shapely import Point, Polygon
from shapely import prepare as prepare_polygon

from .._lib.colour import Colour
from .._lib.geometry import Coord2DLike


class AbstractShape(metaclass=ABCMeta):
    """Abstract Shape Type Class"""

    def __init__(self, xy: Coord2DLike, attribute: Any) -> None:
        if len(xy) != 2:
            raise ValueError(
                "size has be an list of two numerals (width, height)")
        self._xy: Tuple[float, float] = tuple(xy)  # type: ignore
        self._attribute = attribute

    @property
    def xy(self) -> Tuple[float, float]:
        return self._xy

    @property
    def attribute(self) -> Any:
        return self._attribute

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
    def make_polygon(self, buffer: int = 0, prepare: bool = True) -> Polygon:
        """make a polygon"""


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

    def variant(
        self,
        xy: Optional[Coord2DLike] = None,
        diameter: Optional[float] = None,
        attribute: Optional[Any] = None,
    ) -> Dot:
        """return a variant of this Rectangle"""
        if xy is None:
            xy = self._xy
        if diameter is None:
            diameter = self._diameter
        if attribute is None:
            attribute = self._attribute
        return Dot(xy, diameter, attribute)

    def make_polygon(self, buffer: int = 0, prepare: bool = True) -> Polygon:
        """returns  polygon"""
        bf = self._diameter / 2 + buffer
        poly = Point(self._xy).buffer(bf, quad_segs=Dot.QUAD_SEGS)
        if prepare:
            prepare_polygon(poly)  # TODO needed?
        return poly

    def __repr__(self):
        return (
            f"Dot(xy={self._xy}, diameter={self._diameter}, "
            + f"attribute = {self._attribute})"
        )


class AbstractRectangle(AbstractShape):
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
        self._size: Tuple[float, float] = tuple(size)  # type: ignore

    @property
    def size(self) -> Tuple[float, float]:
        return self._size

    def make_polygon(self, buffer: int = 0, prepare: bool = True) -> Polygon:
        l = self._xy[0] - self._size[0] / 2
        r = self._xy[0] + self._size[0] / 2
        t = self._xy[1] + self._size[1] / 2
        b = self._xy[1] - self._size[1] / 2
        poly = Polygon(((l, b), (l, t), (r, t), (r, b)))
        if buffer != 0:
            poly = poly.buffer(buffer, quad_segs=AbstractRectangle.QUAD_SEGS)
        if prepare:
            prepare_polygon(poly)  # TODO needed?
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


class Rectangle(AbstractRectangle):

    def __init__(self,
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

        super().__init__(xy, size, attribute)

    def __repr__(self):
        return (f"Rectangle(xy={self._xy}, size={self._size}, "
                + f"attribute='{self._attribute}')")

    def variant(self,
                xy: Optional[Coord2DLike] = None,
                size: Optional[Coord2DLike] = None,
                attribute: Optional[Any] = None,
                ) -> Rectangle:
        """return a variant of this Rectangle"""
        if xy is None:
            xy = self._xy
        if size is None:
            size = self._size
        if attribute is None:
            attribute = self._attribute
        return Rectangle(xy, size, attribute)


class Picture(AbstractRectangle):
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

    def variant(self,
                xy: Optional[Coord2DLike] = None,
                size: Optional[Coord2DLike] = None,
                filename: Optional[str] = None) -> Picture:
        """return a variant of this Rectangle"""
        if xy is None:
            xy = self._xy
        if size is None:
            size = self._size
        if filename is None:
            filename = self.filename
        return Picture(xy, size, filename)
