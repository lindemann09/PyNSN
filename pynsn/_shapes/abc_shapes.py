
from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from abc import ABCMeta, abstractmethod
from typing import Any, Optional, Union
from typing import Sequence, Tuple

import numpy as np
import shapely
from shapely import affinity
from numpy.typing import NDArray
from shapely import Point, Polygon

from .colour import Colour

INCORRECT_COORDINATE = "xy has be an list of two numerals (x, y)"

Coord2D = Tuple[float, float]
Coord2DLike = Union[Coord2D, Sequence[float], NDArray]


class AbstractPoint(metaclass=ABCMeta):

    def __init__(self, xy: Coord2DLike):
        self._xy = np.asarray(xy)
        if len(self._xy) != 2:
            raise ValueError(INCORRECT_COORDINATE)

    @property
    def xy(self) -> NDArray:
        return self._xy

    @xy.setter
    def xy(self, val: Coord2DLike):
        self._xy = np.asarray(val)
        if len(self._xy) != 2:
            raise ValueError(INCORRECT_COORDINATE)

    @property
    def xy_point(self) -> Point:
        return Point(self._xy)

    @classmethod
    @property
    def name(cls) -> str:
        return str(cls.__name__)

    @abstractmethod
    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        """Distance to another shape

        Note: Returns negative distances only for circular shapes, otherwise
        overlapping distances are 0.
        """

    @abstractmethod
    def dwithin(self, shape: Union[AbstractPoint, AbstractShape], distance: float) -> bool:
        """True if point is in given distance to a shape (dist)

        Using this function is more efficient for non-circular shapes than
        computing the distance and comparing it with dist.
        """

    @abstractmethod
    def is_inside(self, shape: AbstractPoint,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        """True is shapes fully inside the shapes (dist)
        """

    def todict(self) -> dict:
        """dict representation of the object"""
        return {"type": str(self.__class__.__name__),
                "xy": self.xy.tolist()}


class AbstractShape(metaclass=ABCMeta):
    """Abstract Shape Type Class"""

    def __init__(self,
                 xy: Coord2DLike,
                 size: Coord2DLike,
                 attribute: Any) -> None:
        self._xy = np.asarray(xy)
        if len(self._xy) != 2:
            raise ValueError(INCORRECT_COORDINATE)
        self._size = np.asarray(size)
        if len(self._size) != 2:
            raise ValueError(
                "size has be an list of two numerals (width, height)")
        self._attribute = None
        self._polygon = None
        self.attribute = attribute  # setter

    @property
    def xy(self) -> NDArray:
        return self._xy

    @xy.setter
    def xy(self, val: Coord2DLike):
        xy = np.asarray(val)
        if len(xy) != 2:
            raise ValueError(INCORRECT_COORDINATE)

        if isinstance(self._polygon, Polygon):
            move_xy = (xy[0] - self._xy[0], xy[1] - self._xy[1])
            self._polygon = shapely.transform(self._polygon,
                                              lambda x: x + move_xy)
            shapely.prepare(self._polygon)
        self._xy = xy

    @property
    def size(self) -> NDArray:
        return self._size

    @size.setter
    def size(self, val: Coord2DLike):
        """scale the size of the object
        """
        val = np.asarray(val)
        fact = val / self._size
        self._size = val
        if isinstance(self._polygon, Polygon):
            self._polygon = affinity.scale(self._polygon,
                                           xfact=fact[0], yfact=fact[1])
            shapely.prepare(self._polygon)

    @property
    def attribute(self) -> Any:
        return self._attribute

    @attribute.setter
    def attribute(self, val: Any):
        if val is None:
            self._attribute = None
        else:
            try:
                self._attribute = Colour(val)
            except TypeError:
                self._attribute = val

    @property
    def colour(self) -> Colour:
        """colour of the shape
        """
        if isinstance(self._attribute, Colour):
            return self._attribute
        else:
            return Colour(None)

    @property
    def width(self) -> float:
        return self.size[0]

    @property
    def height(self) -> float:
        return self.size[1]

    def scale(self, factor: float):
        """scale the size of the object
        """
        if factor != 1:
            self._size = self._size * factor
            if isinstance(self._polygon, Polygon):
                self._polygon = affinity.scale(self._polygon,
                                               xfact=factor, yfact=factor)
                shapely.prepare(self._polygon)

    def make_polygon(self) -> None:
        """enforce the creation of polygon, if it does not exist yet
        """
        if self._polygon is None:
            self._polygon = self.polygon

    @classmethod
    def name(cls) -> str:
        return str(cls.__name__)

    ## abstract methods ###
    @abstractmethod
    def todict(self) -> dict:
        """dict representation of the object"""
        return {"type": self.name(),
                "xy": self.xy.tolist(),
                "attr": str(self.attribute)}

    @property
    @abstractmethod
    def polygon(self) -> Polygon:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        """"""

    @abstractmethod
    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        """Distance to another shape

        Note: Returns negative distances only for circular shapes, otherwise
        overlapping distances are 0.
        """

    @abstractmethod
    def dwithin(self, shape: Union[AbstractPoint, AbstractShape], distance: float) -> bool:
        """True is shapes are within a given distance

        Using this function is more efficient for non-circular shapes than
        computing the distance and comparing it with dist.
        """

    @abstractmethod
    def is_inside(self, shape: AbstractShape,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        """True is shapes fully inside the shapes (dist)
        """


class AbstractCircularShape(AbstractShape, metaclass=ABCMeta):
    """Abstract Class for Circular Shapes"""
    QUAD_SEGS = 32  # line segments used to approximate dot

    @abstractmethod
    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        pass

    @abstractmethod
    def dwithin(self, shape: AbstractShape, distance: float) -> bool:
        pass

    @abstractmethod
    def is_inside(self, shape: AbstractShape,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        pass


def is_in_shape(a: Union[AbstractPoint, AbstractShape],
                b: AbstractShape,
                b_exterior_ring: Optional[shapely.LinearRing] = None,
                min_dist_boarder: float = 0) -> bool:
    """Returns True if shape or PointType a is fully inside shape b while taking
    into account a minimum to the shape boarder.

    If b_exterior_ring  is not defined, it has to created to determine distance
    to boarder for non-circular shapes. That is, if b_exterior_ring
    is already known in this case, specifying the parameter will improve performance.
    """

    if isinstance(a, AbstractPoint):
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
