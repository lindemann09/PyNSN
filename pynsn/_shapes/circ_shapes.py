"""Point2D and circular shapes

For performance reasons, circular shapes (Dot & Ellipse) are represented merely by
positions and radii. Spatial module between circular shapes not on
shapely.polygons. Polygons will be created if required only.

"""

from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from copy import deepcopy
from typing import Any, Optional, Tuple, Union
from numpy.typing import NDArray

import numpy as np
import shapely
from shapely import Point, Polygon
from shapely.affinity import scale

from ..types import Coord2DLike
from .abc_shapes import CircularShapeType, PointType, ShapeType, is_in_shape


class Point2D(PointType):

    def distance(self, shape: Union[PointType, ShapeType]) -> float:
        """Distance to another shape

        Note: Returns negative distances only for circular shapes, otherwise
        overlapping distances are 0.
        """
        if isinstance(shape, (PointType, CircularShapeType)):
            return _distance_circ_circ(self, shape)

        else:
            return shapely.distance(self.xy_point, shape.polygon)

    def dwithin(self, shape: Union[PointType, ShapeType], distance: float) -> bool:
        """True if point is in given distance to a shape (dist)

        Using this function is more efficient for non-circular shapes than
        computing the distance and comparing it with dist.
        """

        if isinstance(shape, (PointType, CircularShapeType)):
            return self.distance(shape) < distance
        else:
            return shapely.dwithin(self.xy_point, shape.polygon, distance=distance)

    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        """True is shapes fully inside the shapes (dist)
        """
        if isinstance(shape, (Dot, Ellipse)):
            return _is_circ_in_circ(self, b=shape,
                                    min_dist_boarder=min_dist_boarder)
        else:
            return is_in_shape(self,
                               b=shape.polygon,
                               b_exterior_ring=shape_exterior_ring,
                               min_dist_boarder=min_dist_boarder)


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
        self._size = np.asarray(size)
        if len(self._size) != 2:
            raise ValueError(
                "size has be an list of two numerals (width, height)")

    @property
    def size(self) -> NDArray:
        return self._size

    def diameter(self, theta: float) -> float:
        """Returns the diameter at a certain angle (theta)"""
        d = self._size
        return (d[0] * d[1]) / np.sqrt((d[0] * np.sin(theta))**2
                                       + (d[1] * np.cos(theta))**2)

    @property
    def polygon(self) -> Polygon:  # lazy polygon creation
        if self._polygon is None:
            circle = Point(self._xy).buffer(1, quad_segs=Dot.QUAD_SEGS)
            self._polygon = scale(circle, self.size[0]/2, self.size[1]/2)
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

    def distance(self, shape: Union[PointType, ShapeType]) -> float:
        if isinstance(shape, (PointType, CircularShapeType)):
            return _distance_circ_circ(self, shape)

        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: ShapeType, distance: float) -> bool:
        if isinstance(shape, (PointType, CircularShapeType)):
            return self.distance(shape) < distance
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=distance)

    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        if isinstance(shape, CircularShapeType):
            return _is_circ_in_circ(self, b=shape,
                                    min_dist_boarder=min_dist_boarder)
        else:
            return is_in_shape(self, b=shape,
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
    def size(self) -> NDArray:
        return np.array((self._diameter, self._diameter))

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

    def distance(self, shape: Union[PointType, ShapeType]) -> float:
        if isinstance(shape, (PointType, CircularShapeType)):
            return _distance_circ_circ(self, shape)

        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: ShapeType, distance: float) -> bool:
        if isinstance(shape, (PointType, CircularShapeType)):
            return self.distance(shape) < distance
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=distance)

    def is_inside(self, shape: ShapeType,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        if isinstance(shape, CircularShapeType):
            return _is_circ_in_circ(self, b=shape,
                                    min_dist_boarder=min_dist_boarder)
        else:
            return is_in_shape(self, b=shape,
                               b_exterior_ring=shape_exterior_ring,
                               min_dist_boarder=min_dist_boarder)


def _distance_circ_circ(a: Union[PointType, CircularShapeType],
                        b: Union[PointType, CircularShapeType]) -> float:
    """Returns the distance between a circular shape or PointType and other
    circular shape or PointType
    """
    d_xy = np.asarray(a.xy) - b.xy
    theta = None
    if isinstance(a, Dot):
        dia_a = a.diameter
    elif isinstance(a, Ellipse):
        if theta is None:
            theta = np.arctan2(d_xy[1], d_xy[0])
        dia_a = a.diameter(theta=theta)
    elif isinstance(a, PointType):
        dia_a = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(a)}")

    if isinstance(b, Dot):
        dia_b = b.diameter
    elif isinstance(b, Ellipse):
        if theta is None:
            theta = np.arctan2(d_xy[1], d_xy[0])
        dia_b = b.diameter(theta=theta)
    elif isinstance(b, PointType):
        dia_b = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")

    return np.hypot(d_xy[0], d_xy[1]) - (dia_a + dia_b) / 2


def _is_circ_in_circ(a: Union[PointType, CircularShapeType],
                     b: Union[PointType, CircularShapeType],
                     min_dist_boarder: float = 0) -> bool:
    """Returns True if circular shape or PointType is inside another circular
    shape or Point while taking into account a minimum to the shape boarder.
    """
    d_xy = np.asarray(a.xy) - b.xy
    theta = None
    if isinstance(a, Dot):
        a_diameter = a.diameter
    elif isinstance(a, Ellipse):
        if theta is None:
            theta = np.arctan2(d_xy[1], d_xy[0])
        a_diameter = a.diameter(theta=theta)
    elif isinstance(a, PointType):
        a_diameter = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(a)}")

    if isinstance(b, Dot):
        # dot - dot
        b_diameter = b.diameter
    elif isinstance(b, Ellipse):
        # dot/ellipse - dot/ellipse
        if theta is None:
            theta = np.arctan2(d_xy[1], d_xy[0])
        b_diameter = b.diameter(theta=theta)
    elif isinstance(b, PointType):
        b_diameter = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")

    max_ctr_dist = (b_diameter - a_diameter) / 2 - min_dist_boarder
    return max_ctr_dist > 0 and np.hypot(d_xy[0], d_xy[1]) < max_ctr_dist
