
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from typing import Any, Optional, Tuple, Union

import numpy as np
import shapely
from numpy.typing import NDArray
from shapely import Point, Polygon

from ..types import Coord2D, Coord2DLike
from .colour import Colour

INCORRECT_COORDINATE = "xy has be an list of two numerals (x, y)"


class PointType(metaclass=ABCMeta):

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

    @property
    def xy_point(self) -> Point:
        return Point(self._xy)


    @abstractmethod
    def distance(self, shape: Union[PointType, ShapeType]) -> float:
            """Distance to another shape

            Note: Returns negative distances only for circular shapes, otherwise
            overlapping distances are 0.
            """


    @abstractmethod
    def dwithin(self, shape: Union[PointType, ShapeType], dist: float) -> bool:
        """True if point is in given distance to a shape (dist)

        Using this function is more efficient for non-circular shapes than
        computing the distance and comparing it with dist.
        """

    @abstractmethod
    def is_inside(self, shape: PointType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        """True is shapes fully inside the shapes (dist)
        """



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
    def distance(self, shape: Union[PointType, ShapeType]) -> float:
        """Distance to another shape

        Note: Returns negative distances only for circular shapes, otherwise
        overlapping distances are 0.
        """

    @abstractmethod
    def dwithin(self, shape: Union[PointType, ShapeType], dist: float) -> bool:
        """True is shapes are within a given distance (dist)

        Using this function is more efficient for non-circular shapes than
        computing the distance and comparing it with dist.
        """

    @abstractmethod
    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        """True is shapes fully inside the shapes (dist)
        """


class CircularShapeType(ShapeType, metaclass=ABCMeta):
    """Abstract Class for Circular Shapes"""
    QUAD_SEGS = 32  # line segments used to approximate dot

    @abstractmethod
    def distance(self, shape: Union[PointType, ShapeType]) -> float:
        pass

    @abstractmethod
    def dwithin(self, shape: ShapeType, distance: float) -> bool:
        pass

    @abstractmethod
    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        pass

def is_in_shape(a: Union[PointType, ShapeType],
                b: ShapeType,
                b_exterior_ring: Optional[shapely.LinearRing] = None,
                min_dist_boarder: float = 0) -> bool:
    """Returns True if shape or PointType a is fully inside shape b while taking
    into account a minimum to the shape boarder.

    If b_exterior_ring  is not defined, it has to created to determine distance
    to boarder for non-circular shapes. That is, if b_exterior_ring
    is already known in this case, specifying the parameter will improve performance.
    """

    if isinstance(a, PointType):
        a_polygon = a.xy_point
    else:
        a_polygon = a.polygon

    if not shapely.contains_properly(b.polygon, a_polygon):
        # target is not inside if fully covered
        return False
    else:
        if min_dist_boarder > 0:
            # is target to too close to target_area_ring -> False
            if not isinstance(b_exterior_ring, shapely.LinearRing):
                b_exterior_ring = shapely.get_exterior_ring(b.polygon)
            return not shapely.dwithin(a_polygon, b_exterior_ring,
                                       distance=min_dist_boarder)
        else:
            return True
