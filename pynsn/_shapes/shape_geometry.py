"""
calculating spatial relations optimized for circular shapes

Note: The module based the calculation of relations between circular shapes
not on shapely.polygons and used merely positions and radii.
"""

from typing import Optional

import numpy as np
import shapely
from numpy.typing import NDArray

from .shapes import CircularShapeType, Dot, Ellipse, ShapeType, Point2D


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

# FIXME WEITER HIER: move is_inside to shapes


def distance_point_point(a: Point2D, b: Point2D):
    d_xy = np.asarray(a.xy) - b.xy
    return np.hypot(d_xy[0], d_xy[1])


def distance_point_circ(a: Point2D, b: CircularShapeType):
    d_xy = np.asarray(a.xy) - b.xy
    if isinstance(b, Dot):
        dia = b.diameter
    elif isinstance(b, Ellipse):
        theta = np.arctan2(d_xy[0], d_xy[1])
        dia = b.diameter(theta=theta)
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")
    # center dist - radius_b
    return np.hypot(d_xy[0], d_xy[1]) - dia / 2


def distance_circ_circ(a: CircularShapeType, b: CircularShapeType) -> float:
    d_xy = np.asarray(a.xy) - b.xy
    theta = np.arctan2(d_xy[0], d_xy[1])
    if isinstance(a, Dot):
        dia_a = a.diameter
    elif isinstance(a, Ellipse):
        dia_a = a.diameter(theta=theta)
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(a)}")
    if isinstance(b, Dot):
        dia_b = b.diameter
    elif isinstance(b, Ellipse):
        dia_b = b.diameter(theta=theta)
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(b)}")

    return np.hypot(d_xy[0], d_xy[1]) - (dia_a + dia_b) / 2


def is_inside(a: ShapeType,
              b: ShapeType,
              b_exterior_ring: Optional[shapely.LinearRing] = None,
              min_dist_boarder: int = 0) -> bool:
    """Returns True if shape a is fully inside shape b

    If b_exterior_ring  is not defined, it has to created to determine distance
    to boarder for non-circular shapes. That is, if b_exterior_ring
    is already known in this case, specifying the parameter will improve performance.
    """
    if isinstance(a, Dot) and isinstance(b, Dot):
        # dot - dot
        center_dist = np.hypot(b.xy[0] - a.xy[0], b.xy[1] - a.xy[1])
        max_ctr_dist = (b.diameter - a.diameter) / 2 - min_dist_boarder
        return max_ctr_dist > 0 and center_dist < max_ctr_dist
    elif isinstance(a, CircularShapeType) and isinstance(b, CircularShapeType):
        # dot/ellipse - dot/ellipse
        center_dist = np.hypot(b.xy[0] - a.xy[0], b.xy[1] - a.xy[1])
        raise NotImplementedError()
    else:
        # based on shapes
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


def distance_point_dot_array(p: Point2D,
                             dots_xy: NDArray, dots_diameter: NDArray) -> NDArray[np.float_]:
    """Distances of Point2D to multiple dots
    """
    d_xy = dots_xy - p.xy
    return np.hypot(d_xy[:, 0], d_xy[:, 1]) - dots_diameter / 2


def distance_point_ellipses_array(p: Point2D,
                                  ellipses_xy: NDArray, ellipse_sizes: NDArray) -> NDArray[np.float_]:
    """Distances of Point2D to multiple dots
    """
    # coordinate -> ellipses
    d_xy = ellipses_xy - p.xy
    theta = np.arctan2(d_xy[:, 0], d_xy[:, 1])
    # radii of ellipses in array to circ_shape
    ellipse_dia = ellipse_diameter(size=ellipse_sizes, theta=theta)
    # center dist - radius_a - radius_b
    return np.hypot(d_xy[:, 0], d_xy[:, 1]) - ellipse_dia / 2


def distance_circ_dot_array(circ_shape: CircularShapeType,
                            dots_xy: NDArray, dots_diameter: NDArray) -> NDArray[np.float_]:
    """Distances circular shape to multiple dots
    """
    d_xy = dots_xy - circ_shape.xy
    if isinstance(circ_shape, Dot):
        circ_dia = circ_shape.diameter
    elif isinstance(circ_shape, Ellipse):
        circ_dia = ellipse_diameter(
            size=np.atleast_2d(circ_shape.size),
            theta=np.arctan2(d_xy[:, 1], d_xy[:, 0]))  # ellipse radius to each dot in the array
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(circ_shape)}")

    # center dist - radius_a - radius_b
    return np.hypot(d_xy[:, 0], d_xy[:, 1]) - (dots_diameter + circ_dia) / 2


def distance_circ_ellipse_array(circ_shape: CircularShapeType,
                                ellipses_xy: NDArray, ellipse_sizes: NDArray) -> NDArray[np.float_]:
    """Distance circular shape to multiple ellipses
    """
    d_xy = ellipses_xy - circ_shape.xy
    theta = np.arctan2(d_xy[:, 0], d_xy[:, 1])
    # radii of ellipses in array to circ_shape
    ellipse_dia = ellipse_diameter(size=ellipse_sizes, theta=theta)
    if isinstance(circ_shape, Dot):
        shape_dia = circ_shape.diameter
    elif isinstance(circ_shape, Ellipse):
        shape_dia = ellipse_diameter(size=np.atleast_2d(circ_shape.size),
                                     theta=theta)  # ellipse radius to each ellipse in the array
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(circ_shape)}")

    # center dist - radius_a - radius_b
    return np.hypot(d_xy[:, 0], d_xy[:, 1]) - (ellipse_dia + shape_dia)/2
