"""
calculating spatial relations optimized for circular shapes

Note: The module based the calculation of relations between circular shapes
not on shapely.polygons and used merely positions and radii.
"""

from typing import Optional

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


def distance_point_point(a: _shapes.Point2D, b: _shapes.Point2D) -> float:
    return np.hypot(b.xy[0] - a.xy[0], b.xy[1] - a.xy[1])


def distance_point_circ(a: _shapes.Point2D, b: _shapes.CircularShapeType) -> float:
    d_xy = np.asarray(a.xy) - b.xy
    if isinstance(b, _shapes.Dot):
        dia = b.diameter
    elif isinstance(b, _shapes.Ellipse):
        theta = np.arctan2(d_xy[0], d_xy[1])
        dia = b.diameter(theta=theta)
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")
    # center dist - radius_b
    return np.hypot(d_xy[0], d_xy[1]) - dia / 2


def distance_dot_dot(a: _shapes.Dot, b: _shapes.Dot) -> float:
    return np.hypot(b.xy[0] - a.xy[0], b.xy[1] - a.xy[1]) - (a.diameter + b.diameter) / 2


def distance_circ_circ(a: _shapes.CircularShapeType, b: _shapes.CircularShapeType) -> float:
    d_xy = np.asarray(a.xy) - b.xy
    theta = np.arctan2(d_xy[0], d_xy[1])
    if isinstance(a, _shapes.Dot):
        dia_a = a.diameter
    elif isinstance(a, _shapes.Ellipse):
        dia_a = a.diameter(theta=theta)
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(a)}")
    if isinstance(b, _shapes.Dot):
        dia_b = b.diameter
    elif isinstance(b, _shapes.Ellipse):
        dia_b = b.diameter(theta=theta)
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")

    return np.hypot(d_xy[0], d_xy[1]) - (dia_a + dia_b) / 2


def is_point_in_circ(p: _shapes.Point2D,
                     shape: _shapes.CircularShapeType,
                     min_dist_boarder: float = 0) -> bool:
    """Returns True if point a is inside circular shape"""
    d_xy = np.asarray(shape.xy) - p.xy
    if isinstance(shape, _shapes.Dot):
        # dot - dot
        return np.hypot(d_xy[0], d_xy[1]) < shape.diameter/2 - min_dist_boarder
    elif isinstance(shape, _shapes.Ellipse):
        # dot/ellipse - dot/ellipse
        theta = np.arctan2(d_xy[0], d_xy[1])
        return np.hypot(d_xy[0], d_xy[1]) < shape.diameter(theta=theta)/2 - min_dist_boarder
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(shape)}")


def is_dot_in_dot(a: _shapes.Dot,
                  b: _shapes.Dot,
                  min_dist_boarder: float = 0) -> bool:
    """Returns True if point a is inside circular shape"""
    d_xy = np.asarray(b.xy) - a.xy
    max_ctr_dist = (b.diameter - a .diameter) / 2 - min_dist_boarder
    return max_ctr_dist > 0 and np.hypot(d_xy[0], d_xy[1]) < max_ctr_dist


def is_circ_in_circ(a: _shapes.CircularShapeType,
                    b: _shapes.CircularShapeType,
                    min_dist_boarder: float = 0) -> bool:
    """Returns True if point a is inside circular shape"""
    d_xy = np.asarray(a.xy) - b.xy
    theta = np.arctan2(d_xy[0], d_xy[1])
    if isinstance(a, _shapes.Dot):
        a_diameter = a.diameter
    elif isinstance(a, _shapes.Ellipse):
        a_diameter = a.diameter(theta=theta)
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(a)}")

    if isinstance(b, _shapes.Dot):
        # dot - dot
        max_ctr_dist = (b.diameter - a_diameter) / 2 - min_dist_boarder
    elif isinstance(b, _shapes.Ellipse):
        # dot/ellipse - dot/ellipse
        max_ctr_dist = (b.diameter(theta=theta) -
                        a_diameter) / 2 - min_dist_boarder
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")

    return max_ctr_dist > 0 and np.hypot(d_xy[0], d_xy[1]) < max_ctr_dist


def is_point_in_shape(a: _shapes.Point2D,
                      b: _shapes.ShapeType,
                      b_exterior_ring: Optional[shapely.LinearRing] = None,
                      min_dist_boarder: float = 0) -> bool:
    """Returns True if point a is fully inside shape b

    If b_exterior_ring  is not defined, it has to created to determine distance
    to boarder for non-circular shapes. That is, if b_exterior_ring
    is already known in this case, specifying the parameter will improve performance.
    """

    if not shapely.contains_properly(b.polygon, a.xy_point):
        # target is not inside if fully covered
        return False
    else:
        if min_dist_boarder > 0:
            # is target to too close to target_area_ring -> False
            if not isinstance(b_exterior_ring, shapely.LinearRing):
                b_exterior_ring = shapely.get_exterior_ring(b.polygon)
            return not shapely.dwithin(a.xy_point, b_exterior_ring,
                                       distance=min_dist_boarder)
        else:
            return True


def is_shape_in_shape(a: _shapes.ShapeType,
                      b: _shapes.ShapeType,
                      b_exterior_ring: Optional[shapely.LinearRing] = None,
                      min_dist_boarder: float = 0) -> bool:
    """Returns True if shape a is fully inside shape b

    If b_exterior_ring  is not defined, it has to created to determine distance
    to boarder for non-circular shapes. That is, if b_exterior_ring
    is already known in this case, specifying the parameter will improve performance.
    """

    if not shapely.contains_properly(b.polygon, a.polygon):
        # target is not inside if fully covered
        return False
    else:
        if min_dist_boarder > 0:
            # is target to too close to target_area_ring -> False
            if not isinstance(b_exterior_ring, shapely.LinearRing):
                b_exterior_ring = shapely.get_exterior_ring(b.polygon)
            return not shapely.dwithin(a.polygon, b_exterior_ring,
                                       distance=min_dist_boarder)
        else:
            return True
