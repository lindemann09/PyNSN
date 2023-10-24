"""Shape Classes

Rectanglur shapes and Polygons are based on shapely.Polygon.
For performance reasons, circular shapes (Dot & Ellipse) are represented merely by
positions and radii. Polygons will be created if required only.

Note: For performance reasons, spatial relations module based the calculation of
relations between circular shapes not on shapely.polygons.
"""

from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from os import path
from typing import Any, Tuple, Union
from copy import deepcopy

import shapely
from shapely.affinity import scale
from shapely import Polygon, Point

from .._lib.colour import Colour
from .._lib.geometry import Coord2DLike, Coord2D


class ShapeType(metaclass=ABCMeta):
    """Abstract Shape Type Class"""
    ID = -1

    def __init__(self,
                 xy: Coord2DLike,
                 attribute: Any) -> None:
        if len(xy) != 2:
            raise ValueError("xy has be an list of two numerals (x, y)")

        self._xy = tuple(xy)
        self._attribute = attribute
        self._polygon = None

    @property
    def xy(self) -> Coord2D:
        return self._xy  # type: ignore

    @xy.setter
    def xy(self, val: Coord2DLike):
        if len(val) != 2:
            raise ValueError("xy has be an list of two numerals (x, y)")

        xy = tuple(val)
        if isinstance(self._polygon, Polygon):
            move_xy = (xy[0] - self._xy[0], xy[1] - self._xy[1])
            self._polygon = shapely.transform(self._polygon,
                                              lambda x: x + move_xy)  # FIXME Are transformed polygons still prepared?
        self._xy = xy

    @property
    def attribute(self) -> Any:
        return self._attribute

    @attribute.setter
    def attribute(self, val: Any):
        self._attribute = val

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

    @property
    def width(self) -> float:
        return self.size[0]

    @property
    def height(self) -> float:
        return self.size[1]

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

    def copy(self):
        """returns a copy of the object
        """
        return deepcopy(self)

    @property
    @abstractmethod
    def size(self) -> Tuple[float, float]:
        pass

    @property
    @abstractmethod
    def polygon(self) -> Polygon:
        pass

    # def delete_polygon(self):
    #     self._polygon = None

    # def set_polygon(self, polygon: Polygon) -> None:
    #     """This method should usually not been used, unless you know what you do.

    #     Polygon is automatically created, if required. Please use the property
    #     `polygon.
    #     """
    #     self._polygon = polygon

    @abstractmethod
    def __repr__(self) -> str:
        """"""


class PolygonShape(ShapeType):
    ID = 5

    def __init__(self, polygon: Polygon, attribute: Any = None):
        ctr = polygon.centroid
        super().__init__(xy=(ctr.x, ctr.y), attribute=attribute)
        shapely.prepare(polygon)
        self._polygon = polygon

    @property
    def polygon(self) -> Polygon:
        return self._polygon

    @property
    def size(self) -> Tuple[float, float]:
        b = shapely.bounds(self._polygon)
        return b[2:4] - b[0:2]  # [bound width, bound height]

    def __repr__(self):
        return (f"PolygonShape(xy={self._xy}, size={self.size}, "
                + f"attribute='{self._attribute}')")


class Ellipse(ShapeType):
    QUAD_SEGS = 32  # line segments used to approximate dot
    ID = 3

    def __init__(self,
                 xy: Coord2DLike,
                 size: Coord2DLike,
                 attribute: Any = None
                 ) -> None:
        """Initialize a dot

        Parameters
        ----------
        xy : tuple of two numeric
        size : x and y diameter
        attribute : attribute (string, optional)
        """
        super().__init__(xy=xy, attribute=attribute)
        if len(size) != 2:
            raise ValueError(
                "size has be an list of two numerals (width, height)")
        self._size = tuple(size)

    @property
    def size(self) -> Coord2D:
        return self._size  # type: ignore

    @property
    def polygon(self) -> Polygon:  # lazy polygon creation
        if self._polygon is None:
            circle = Point(self._xy).buffer(1, quad_segs=Dot.QUAD_SEGS)
            self._polygon = scale(circle, 15, 20)
            shapely.prepare(self._polygon)  # FIXME needed?
        return self._polygon

    def __repr__(self):
        return (
            f"Ellipse(xy={self._xy}, size={self.size}, "
            + f"attribute = {self._attribute})"
        )


class Dot(Ellipse):
    ID = 1

    def __init__(self,
                 xy: Coord2DLike,
                 diameter: float,
                 attribute: Any = None
                 ) -> None:
        """Initialize a dot

        Parameters
        ----------
        xy : tuple of two numeric
        diameter : numeric
        attribute : attribute (string, optional)
        """
        super().__init__(xy=xy,
                         size=(diameter, diameter),
                         attribute=attribute)

    @property
    def diameter(self) -> float:
        return self._size[0]

    @property
    def polygon(self) -> Polygon:  # lazy polygon creation
        if self._polygon is None:
            r = self.diameter / 2
            self._polygon = Point(self._xy).buffer(r, quad_segs=Dot.QUAD_SEGS)
            shapely.prepare(self._polygon)  # FIXME needed?
        return self._polygon

    def __repr__(self):
        return (
            f"Dot(xy={self._xy}, diameter={self.diameter}, "
            + f"attribute = {self._attribute})"
        )


class Rectangle(ShapeType):
    ID = 2

    def __init__(self,
                 xy: Coord2DLike,
                 size: Coord2DLike,
                 attribute: Any = None
                 ):
        super().__init__(xy=xy, attribute=attribute)
        # make polygon
        l = xy[0] - size[0] / 2
        r = xy[0] + size[0] / 2
        t = xy[1] + size[1] / 2
        b = xy[1] - size[1] / 2
        self._polygon = Polygon(((l, b), (l, t), (r, t), (r, b)))
        shapely.prepare(self._polygon)  # FIXME needed?

    @property
    def polygon(self) -> Polygon:
        return self._polygon

    @property
    def size(self):
        lb = shapely.get_coordinates(self._polygon)[0]  # left bottom
        rt = shapely.get_coordinates(self._polygon)[2]  # right top
        return rt - lb

    @property
    def proportion(self) -> float:
        """Proportion of the rectangle (width/height)"""
        return self.size[0] / self.size[1]

    @property
    def left_bottom(self) -> Coord2D:
        return shapely.get_coordinates(self._polygon)[0]

    @property
    def left_top(self) -> Coord2D:
        return shapely.get_coordinates(self._polygon)[1]

    @property
    def right_top(self) -> Coord2D:
        return shapely.get_coordinates(self._polygon)[2]

    @property
    def right_bottom(self) -> Coord2D:
        return shapely.get_coordinates(self._polygon)[3]

    def __repr__(self):
        return (f"Rectangle(xy={self._xy}, size={self.size}, "
                + f"attribute='{self._attribute}')")


class Picture(Rectangle):
    ID = 4
    ATTR_PREFIX = "p:"

    def __init__(self, xy: Coord2DLike, size: Coord2DLike, filename: str) -> None:
        """Initialize a Picture

        Rectangle can also consist of a picture

        Parameters
        ----------
        xy : tuple
            tuple of two numeric
        size : tuple
            tuple of two numeric
        filename : str
            path the picture file
        """

        super().__init__(xy=xy, size=size, attribute=Picture.ATTR_PREFIX + filename)

    def __repr__(self):
        return (f"Picture(xy={self._xy}, size={self.size}, " +
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


SHAPE_LABEL = {
    Dot.ID: "Dot",
    Rectangle.ID: "Rectangle",
    Picture.ID: "Picture",
    Ellipse.ID: "Ellipse",
    PolygonShape.ID: "Polygon"}
