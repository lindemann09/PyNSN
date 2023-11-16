"""Point2D and circular shapes

For performance reasons, circular shapes (Dot & Ellipse) are represented merely by
positions and radii. Spatial module between circular shapes not on
shapely.polygons. Polygons will be created if required only.

"""

from __future__ import annotations

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from typing import Any, Optional,  Union
from numpy.typing import NDArray

import numpy as np
import shapely
from shapely import Point, Polygon, affinity

from .abc_shapes import AbstractCircularShape, AbstractPoint, \
    AbstractShape, is_in_shape, Coord2DLike


class Point2D(AbstractPoint):

    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        """Distance to another shape

        Note: Returns negative distances only for circular shapes, otherwise
        overlapping distances are 0.
        """
        if isinstance(shape, (AbstractPoint, AbstractCircularShape)):
            return _distance_circ_circ(self, shape)

        else:
            return shapely.distance(self.xy_point, shape.polygon)

    def dwithin(self, shape: Union[AbstractPoint, AbstractShape], distance: float) -> bool:
        """True if point is in given distance to a shape (dist)

        Using this function is more efficient for non-circular shapes than
        computing the distance and comparing it with dist.
        """

        if isinstance(shape, (AbstractPoint, AbstractCircularShape)):
            return self.distance(shape) < distance
        else:
            return shapely.dwithin(self.xy_point, shape.polygon, distance=distance)

    def is_inside(self, shape: AbstractShape,
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


class Ellipse(AbstractCircularShape):

    def __init__(self,
                 size: Coord2DLike,
                 xy: Coord2DLike = (0, 0),
                 attribute: Any = None
                 ) -> None:
        """Initialize a dot

        Parameters
        ----------
        xy : tuple of two numeric
        size : x and y diameter
        attribute : attribute (string, optional)
        """
        super().__init__(size=size, xy=xy,  attribute=attribute)

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
            self._polygon = affinity.scale(
                circle, self.size[0]/2, self.size[1]/2)
            shapely.prepare(self._polygon)
        return self._polygon

    def __repr__(self):
        return (
            f"Ellipse(xy={self.xy}, size={self.size}, "
            + f"attribute = {self.attribute})"
        )

    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        if isinstance(shape, (AbstractPoint, AbstractCircularShape)):
            return _distance_circ_circ(self, shape)

        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: AbstractShape, distance: float) -> bool:
        if isinstance(shape, (AbstractPoint, AbstractCircularShape)):
            return self.distance(shape) < distance
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=distance)

    def is_inside(self, shape: AbstractShape,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        if isinstance(shape, AbstractCircularShape):
            return _is_circ_in_circ(self, b=shape,
                                    min_dist_boarder=min_dist_boarder)
        else:
            return is_in_shape(self, b=shape,
                               b_exterior_ring=shape_exterior_ring,
                               min_dist_boarder=min_dist_boarder)

    def todict(self) -> dict:
        d = super().todict()
        d.update({"size": self.size.tolist()})
        return d


class Dot(AbstractCircularShape):

    def __init__(self,
                 diameter: float,
                 xy: Coord2DLike = (0, 0),
                 attribute: Any = None
                 ) -> None:
        """Initialize a dot

        Parameters
        ----------
        xy : tuple of two numeric
        diameter : numeric
        attribute : attribute (string, optional)
        """
        super().__init__(size=np.array((diameter, diameter)), xy=xy,
                         attribute=attribute)

    @property
    def diameter(self) -> float:
        return self._size[0]

    @property
    def polygon(self) -> Polygon:  # lazy polygon creation
        if self._polygon is None:
            r = self._size[0] / 2
            self._polygon = Point(self._xy).buffer(r, quad_segs=Dot.QUAD_SEGS)
            shapely.prepare(self._polygon)  # TODO needed?
        return self._polygon

    def __repr__(self):
        return (
            f"Dot(xy={self._xy}, diameter={self.diameter}, "
            + f"attribute = {self._attribute})"
        )

    def distance(self, shape: Union[AbstractPoint, AbstractShape]) -> float:
        if isinstance(shape, (AbstractPoint, AbstractCircularShape)):
            return _distance_circ_circ(self, shape)

        else:
            return shapely.distance(self.polygon, shape.polygon)

    def dwithin(self, shape: AbstractShape, distance: float) -> bool:
        if isinstance(shape, (AbstractPoint, AbstractCircularShape)):
            return self.distance(shape) < distance
        else:
            return shapely.dwithin(self.polygon, shape.polygon, distance=distance)

    def is_inside(self, shape: AbstractShape,
                  shape_exterior_ring: Optional[shapely.LinearRing] = None,
                  min_dist_boarder: float = 0) -> bool:
        if isinstance(shape, AbstractCircularShape):
            return _is_circ_in_circ(self, b=shape,
                                    min_dist_boarder=min_dist_boarder)
        else:
            return is_in_shape(self, b=shape,
                               b_exterior_ring=shape_exterior_ring,
                               min_dist_boarder=min_dist_boarder)

    def todict(self) -> dict:
        d = super().todict()
        d.update({"diameter": self.diameter})
        return d


def _distance_circ_circ(a: Union[AbstractPoint, AbstractCircularShape],
                        b: Union[AbstractPoint, AbstractCircularShape]) -> float:
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
    elif isinstance(a, AbstractPoint):
        dia_a = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(a)}")

    if isinstance(b, Dot):
        dia_b = b.diameter
    elif isinstance(b, Ellipse):
        if theta is None:
            theta = np.arctan2(d_xy[1], d_xy[0])
        dia_b = b.diameter(theta=theta)
    elif isinstance(b, AbstractPoint):
        dia_b = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")

    return np.hypot(d_xy[0], d_xy[1]) - (dia_a + dia_b) / 2


def _is_circ_in_circ(a: Union[AbstractPoint, AbstractCircularShape],
                     b: Union[AbstractPoint, AbstractCircularShape],
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
    elif isinstance(a, AbstractPoint):
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
    elif isinstance(b, AbstractPoint):
        b_diameter = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")

    max_ctr_dist = (b_diameter - a_diameter) / 2 - min_dist_boarder
    return max_ctr_dist > 0 and np.hypot(d_xy[0], d_xy[1]) < max_ctr_dist
