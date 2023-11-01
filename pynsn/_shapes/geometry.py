"""
calculating spatial relations optimized for circular shapes

Note: The module based the calculation of relations between circular shapes
not on shapely.polygons and used merely positions and radii.
"""

from typing import Optional, Union

import numpy as np
from numpy.typing import NDArray
import shapely

from .. import _shapes


def ellipse_perimeter(sizes: NDArray) -> NDArray[np.float_]:
    """Ramanujan's second approximation of the ellipse perimeter"""
    s = np.atleast_2d(sizes)
    a = s[:, 0]
    b = s[:, 1]
    return np.pi * ((a+b) + (3*(a-b)**2) / (10*(a+b) + np.sqrt(a**2 + 14*a*b + b**2)))


def ellipse_diameter(size: NDArray, theta: NDArray) -> NDArray[np.float_]:
    """Ellipse diameter at a certain angle

    Parameter
    ---------
        Size: NDarray
            2d array with (semi-majors, semi-minors)
        theta: float
            angle in radians
    """
    d = np.atleast_2d(size)
    return (d[:, 0] * d[:, 1]) / np.sqrt((d[:, 0] * np.sin(theta))**2
                                         + (d[:, 1] * np.cos(theta))**2)


def distance_circ_circ(a: Union[_shapes.Point2D, _shapes.CircularShapeType],
                       b: Union[_shapes.Point2D, _shapes.CircularShapeType]) -> float:
    """Returns the distance between a circular shape or Point2D and other
    circular shape or Point2D
    """
    d_xy = np.asarray(a.xy) - b.xy
    theta = None
    if isinstance(a, _shapes.Dot):
        dia_a = a.diameter
    elif isinstance(a, _shapes.Ellipse):
        if theta is None:
            theta = np.arctan2(d_xy[0], d_xy[1])
        dia_a = a.diameter(theta=theta)
    elif isinstance(a, _shapes.Point2D):
        dia_a = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(a)}")

    if isinstance(b, _shapes.Dot):
        dia_b = b.diameter
    elif isinstance(b, _shapes.Ellipse):
        if theta is None:
            theta = np.arctan2(d_xy[0], d_xy[1])
        dia_b = b.diameter(theta=theta)
    elif isinstance(b, _shapes.Point2D):
        dia_b = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")

    return np.hypot(d_xy[0], d_xy[1]) - (dia_a + dia_b) / 2


def is_circ_in_circ(a: Union[_shapes.Point2D, _shapes.CircularShapeType],
                    b: Union[_shapes.Point2D, _shapes.CircularShapeType],
                    min_dist_boarder: float = 0) -> bool:
    """Returns True if circular shape or Point2d is inside another circular
    shape or Point while taking into account a minimum to the shape boarder.
    """
    d_xy = np.asarray(a.xy) - b.xy
    theta = None
    if isinstance(a, _shapes.Dot):
        a_diameter = a.diameter
    elif isinstance(a, _shapes.Ellipse):
        if theta is None:
            theta = np.arctan2(d_xy[0], d_xy[1])
        a_diameter = a.diameter(theta=theta)
    elif isinstance(a, _shapes.Point2D):
        a_diameter = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(a)}")

    if isinstance(b, _shapes.Dot):
        # dot - dot
        b_diameter = b.diameter
    elif isinstance(b, _shapes.Ellipse):
        # dot/ellipse - dot/ellipse
        if theta is None:
            theta = np.arctan2(d_xy[0], d_xy[1])
        b_diameter = b.diameter(theta=theta)
    elif isinstance(b, _shapes.Point2D):
        b_diameter = 0
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")

    max_ctr_dist = (b_diameter - a_diameter) / 2 - min_dist_boarder
    return max_ctr_dist > 0 and np.hypot(d_xy[0], d_xy[1]) < max_ctr_dist


def is_in_shape(a: Union[_shapes.Point2D, _shapes.ShapeType],
                b: _shapes.ShapeType,
                b_exterior_ring: Optional[shapely.LinearRing] = None,
                min_dist_boarder: float = 0) -> bool:
    """Returns True if shape or Point2D a is fully inside shape b while taking
    into account a minimum to the shape boarder.

    If b_exterior_ring  is not defined, it has to created to determine distance
    to boarder for non-circular shapes. That is, if b_exterior_ring
    is already known in this case, specifying the parameter will improve performance.
    """

    if isinstance(a, _shapes.Point2D):
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
