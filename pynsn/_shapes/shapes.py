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

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import Point, Polygon
from shapely.affinity import scale
from .. import _shapes
from ..types import Coord2D, Coord2DLike
from .colour import Colour

INCORRECT_COORDINATE = "xy has be an list of two numerals (x, y)"


class Point2D(object):

    def __init__(self, xy: Coord2DLike):
        if len(xy) != 2:
            raise ValueError(INCORRECT_COORDINATE)
        self._xy = np.asarray(xy)

    @property
    def xy(self) -> NDArray:
        return self._xy

    @xy.setter
    def xy(self, val: Coord2DLike):
        if len(val) != 2:
            raise ValueError(INCORRECT_COORDINATE)
        self._xy = np.asarray(val)

    def xy_point(self) -> Point:
        return Point(self._xy)

    def distance(self, shape: Union[Point2D, ShapeType]) -> float:
        """Distance to another shape

        Note: Returns negative distances only for circular shapes, otherwise
        overlapping distances are 0.
        """
        if isinstance(shape, Point2D):
            return _shapes.shape_geometry.distance_point_point(self, shape)

        elif isinstance(shape, CircularShapeType):
            return _shapes.shape_geometry.distance_point_circ(self, shape)

        else:
            return shapely.distance(self.xy_point, shape.polygon)

    def dwithin(self, shape: Union[Point2D, ShapeType], dist: float) -> bool:
        """True if point is in given distance to a shape (dist)

        Using this function is more efficient for non-circular shapes than
        computing the distance and comparing it with dist.
        """

        if isinstance(shape, Point2D):
            return shapely.dwithin(self.xy_point, shape.xy_point, distance=dist)
        else:
            return shapely.dwithin(self.xy_point, shape.polygon, distance=dist)

    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> float:
        """True is shapes fully inside the shapes (dist)
        """
        if isinstance(shape, (Dot, Ellipse)):
            return _shapes.shape_geometry.is_point_in_circ(self, shape=shape,
                                                           min_dist_boarder=min_dist_boarder)
        else:
            return _shapes.shape_geometry.is_point_in_shape(self,
                                                            b=shape.polygon,
                                                            b_exterior_ring=shape_exterior_ring,
                                                            min_dist_boarder=min_dist_boarder)


class ShapeType(metaclass=ABCMeta):
    """Abstract Shape Type Class"""
    ID = -1

    def __init__(self,
                 xy: Coord2DLike,
                 attribute: Any) -> None:
        if len(xy) != 2:
            raise ValueError(INCORRECT_COORDINATE)
        self._xy = tuple(xy)
        self._attribute = attribute
        self._polygon = None

    @property
    def xy(self) -> Coord2D:
        return self._xy  # type: ignore

    @xy.setter
    def xy(self, val: Coord2DLike):
        if len(val) != 2:
            raise ValueError(INCORRECT_COORDINATE)

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

    def make_polygon(self) -> None:
        """enforce the creation of polygon, if it does not exist yet
        """
        if self._polygon is None:
            self._polygon = self.polygon

    ## abstract methods ###

    @abstractmethod
    def copy(self, new_xy: Optional[Coord2DLike] = None,
             copy_polygon: bool = True) -> ShapeType:
        """returns a copy of the object
        """

    @property
    @abstractmethod
    def size(self) -> Tuple[float, float]:
        pass

    @property
    @abstractmethod
    def polygon(self) -> Polygon:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        """"""

    @abstractmethod
    def distance(self, shape: Union[Point2D, ShapeType]) -> float:
        """Distance to another shape

        Note: Returns negative distances only for circular shapes, otherwise
        overlapping distances are 0.
        """

    @abstractmethod
    def dwithin(self, shape: Union[Point2D, ShapeType], dist: float) -> bool:
        """True is shapes are within a given distance (dist)

        Using this function is more efficient for non-circular shapes than
        computing the distance and comparing it with dist.
        """

    @abstractmethod
    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> float:
        """True is shapes fully inside the shapes (dist)
        """


class CircularShapeType(ShapeType, metaclass=ABCMeta):
    """Abstract Class for Circular Shapes"""
    QUAD_SEGS = 32  # line segments used to approximate dot

    def distance(self, shape: Union[Point2D, ShapeType]) -> float:
        if isinstance(shape, Point2D):
            return _shapes.shape_geometry.distance_point_circ(shape, self)

        elif isinstance(shape, CircularShapeType):
            return _shapes.shape_geometry.distance_circ_circ(self, shape)

        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: ShapeType, dist: float) -> bool:
        if isinstance(shape, (Point2D, CircularShapeType)):
            return self.distance(shape) < dist
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=dist)


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

    def diameter(self, theta: float) -> float:
        """Returns the diameter at a certain angle (theta)"""
        d = self._size
        return (d[0] * d[1]) / np.sqrt((d[0] * np.sin(theta))**2
                                       + (d[1] * np.cos(theta))**2)

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

    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> float:
        if isinstance(shape, (Dot, Ellipse)):
            return _shapes.shape_geometry.is_circ_in_circ(self, b=shape,
                                                          min_dist_boarder=min_dist_boarder)
        else:
            return _shapes.shape_geometry.is_shape_in_shape(self, b=shape,
                                                            b_exterior_ring=shape_exterior_ring,
                                                            min_dist_boarder=min_dist_boarder)


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

    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> float:
        if isinstance(shape, Dot):
            return _shapes.shape_geometry.is_dot_in_dot(self, b=shape,
                                                        min_dist_boarder=min_dist_boarder)
        elif isinstance(shape, Dot):
            return _shapes.shape_geometry.is_circ_in_circ(self, b=shape,
                                                          min_dist_boarder=min_dist_boarder)
        else:
            return _shapes.shape_geometry.is_shape_in_shape(self, b=shape,
                                                            b_exterior_ring=shape_exterior_ring,
                                                            min_dist_boarder=min_dist_boarder)


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

    def distance(self, shape: Union[Point2D, ShapeType]) -> float:
        if isinstance(shape, Point2D):
            return shapely.distance(self.polygon, shape.xy_point)
        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: Union[Point2D, ShapeType], dist: float) -> bool:
        if isinstance(shape, Point2D):
            return shapely.dwithin(self.polygon, shape.xy_point, distance=dist)
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=dist)

    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> float:
        return _shapes.shape_geometry.is_shape_in_shape(self, b=shape,
                                                        b_exterior_ring=shape_exterior_ring,
                                                        min_dist_boarder=min_dist_boarder)


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

    def distance(self, shape: Union[Point2D, ShapeType]) -> float:
        if isinstance(shape, Point2D):
            return shapely.distance(self.polygon, shape.xy_point)
        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: Union[Point2D, ShapeType], dist: float) -> bool:
        if isinstance(shape, Point2D):
            return shapely.dwithin(self.polygon, shape.xy_point, distance=dist)
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=dist)

    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> float:
        return _shapes.shape_geometry.is_shape_in_shape(self, b=shape,
                                                        b_exterior_ring=shape_exterior_ring,
                                                        min_dist_boarder=min_dist_boarder)
