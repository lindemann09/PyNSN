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
from copy import deepcopy
from pathlib import Path
from typing import Any, Optional, Tuple, Union

import shapely
import numpy as np
from numpy.typing import NDArray
from shapely import Point, Polygon
from shapely.affinity import scale

from .._lib.colour import Colour
from .._lib.geometry import Coord2D, Coord2DLike


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
                                              lambda x: x + move_xy)
            shapely.prepare(self._polygon)
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

    @abstractmethod
    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> ShapeType:
        """returns a copy of the object
        """

    def make_polygon(self) -> None:
        """enforce the creation of polygon, if it does not exist yet
        """
        if self._polygon is None:
            self._polygon = self.polygon

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


class CircularShapeType(ShapeType):
    QUAD_SEGS = 32  # line segments used to approximate dot


class Ellipse(CircularShapeType):
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
            f"Ellipse(xy={self.xy}, size={self.size}, "
            + f"attribute = {self.attribute})"
        )

    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> Ellipse:

        if copy_polygon:
            new = deepcopy(self)
            if new_xy is not None:
                new.xy = new_xy
            return new
        elif new_xy is None:
            return Ellipse(xy=self._xy, size=self._size, attribute=self._attribute)
        else:
            return Ellipse(xy=new_xy, size=self._size, attribute=self._attribute)


class Dot(CircularShapeType):
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
                         attribute=attribute)
        self._diameter = diameter

    @property
    def size(self) -> Tuple[float, float]:
        return (self._diameter, self._diameter)

    @property
    def diameter(self) -> float:
        return self._diameter

    @property
    def polygon(self) -> Polygon:  # lazy polygon creation
        if self._polygon is None:
            r = self._diameter / 2
            self._polygon = Point(self._xy).buffer(r, quad_segs=Dot.QUAD_SEGS)
            shapely.prepare(self._polygon)  # FIXME needed?
        return self._polygon

    def __repr__(self):
        return (
            f"Dot(xy={self._xy}, diameter={self._diameter}, "
            + f"attribute = {self._attribute})"
        )

    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> Dot:

        if copy_polygon:
            new = deepcopy(self)
            if new_xy is not None:
                new.xy = new_xy
            return new
        elif new_xy is None:
            return Dot(xy=self._xy, diameter=self.diameter, attribute=self._attribute)
        else:
            return Dot(xy=new_xy, diameter=self.diameter, attribute=self._attribute)


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
    def left_bottom(self) -> NDArray:
        """Returns (left, bottom) as ndarray (x,y)"""
        return shapely.get_coordinates(self._polygon)[0]

    @property
    def right_top(self) -> NDArray:
        """Returns (right, top) as ndarray (x,y)"""
        return shapely.get_coordinates(self._polygon)[2]

    @property
    def left_top(self) -> NDArray[np.float_]:
        """Returns (left, top) as ndarray (x,y)"""
        return shapely.get_coordinates(self._polygon)[1]

    @property
    def right_bottom(self) -> NDArray[np.float_]:
        """Returns (right, bottom) as ndarray (x,y)"""
        return shapely.get_coordinates(self._polygon)[3]

    @property
    def box(self) -> NDArray[np.float_]:
        """Returns (left, bottom, right, top) as NDArray (x0, y0, x1, y1)"""
        return np.append(
            shapely.get_coordinates(self._polygon)[0],
            shapely.get_coordinates(self._polygon)[2])

    def __repr__(self):
        return (f"Rectangle(xy={self._xy}, size={self.size}, "
                + f"attribute='{self._attribute}')")

    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> Rectangle:

        if copy_polygon:
            new = deepcopy(self)
            if new_xy is not None:
                new.xy = new_xy
            return new
        elif new_xy is None:
            return Rectangle(xy=self._xy, size=self.size, attribute=self._attribute)
        else:
            return Rectangle(xy=new_xy, size=self.size, attribute=self._attribute)


class Picture(Rectangle):
    ID = 4

    def __init__(self, xy: Coord2DLike, size: Coord2DLike,
                 path: Union[Path, str]) -> None:
        """Initialize a Picture

        Rectangle can also consist of a picture

        Parameters
        ----------
        xy : tuple
            tuple of two numeric
        size : tuple
            tuple of two numeric
        path : pathlib.path or str
            path the picture file
        """

        super().__init__(xy=xy, size=size, attribute=Path(path))

    def __repr__(self):
        return (f"Picture(xy={self._xy}, size={self.size}, " +
                f"path='{str(self.path)}')")

    @property
    def path(self) -> Path:
        return self._attribute

    def file_exists(self) -> bool:
        """Checks if the file exists"""
        return self.path.isfile()  # type: ignore

    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> Picture:

        if copy_polygon:
            new = deepcopy(self)
            if new_xy is not None:
                new.xy = new_xy
            return new
        elif new_xy is None:
            return Picture(xy=self._xy, size=self.size,
                           path=self.path)
        else:
            return Picture(xy=new_xy, size=self.size,
                           path=self.path)


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

    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> PolygonShape:

        if not copy_polygon:
            raise RuntimeError(
                "copy_polygon = False not is possible for PolygonShape")

        new = deepcopy(self)
        if new_xy is not None:
            new.xy = new_xy
        return new
