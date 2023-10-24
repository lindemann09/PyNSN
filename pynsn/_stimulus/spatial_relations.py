"""
calculating spatial relations optimized for circular shapes

Note: The module based the calculation of relations between circular shapes
not on shapely.polygons and used merely positions and radii.
"""


from numpy.typing import NDArray
import numpy as np
import shapely
from .shape_array import ShapeArray
from .shapes import Dot, Ellipse, Rectangle, Picture, ShapeType


def distance(a: ShapeType, b: ShapeType) -> float:
    if isinstance(a, Dot) and isinstance(b, Dot):
        # dot - dot
        return np.hypot(b.xy[0] - a.xy[0], b.xy[1] - a.xy[1]) \
            - (a.diameter + b.diameter)/2  # -> - radius_a - radius_b
    elif isinstance(a, (Dot, Ellipse)) and isinstance(b, (Dot, Ellipse)):
        # dot/ellipse - dot/ellipse
        raise NotImplemented
    else:
        return shapely.distance(a.polygon, b.polygon)


def dwithin(a: ShapeType, b: ShapeType, dist: float) -> bool:
    if isinstance(a, Dot) and isinstance(b, Dot):
        # dot - dot
        return distance(a, b) <= dist
    elif isinstance(a, (Dot, Ellipse)) and isinstance(b, (Dot, Ellipse)):
        # dot/ellipse - dot/ellipse
        raise NotImplemented
    else:
        return shapely.dwithin(a.polygon, b.polygon, distance=dist)


def contains_properly(a: ShapeType, b: ShapeType) -> bool:  # is inside
    if isinstance(a, Dot) and isinstance(b, Dot):
        # dot - dot
        return distance(a, b) < -2 * b.diameter
    elif isinstance(a, (Dot, Ellipse)) and isinstance(b, (Dot, Ellipse)):
        # dot/ellipse - dot/ellipse
        raise NotImplemented
    else:
        return shapely.contains_properly(a.polygon, b.polygon)


def distance_array(shape: ShapeType, arr: ShapeArray) -> NDArray[float]:
    shapely.prepare(shape.polygon)

    dia = arr.dot_diameter[arr.dot_ids]
    xy = arr.xy[arr.dot_ids, :]


def dwithin_array(shape: ShapeType, arr: ShapeArray, dist: float) -> NDArray[bool]:
    shapely.prepare(shape.polygon)
    pass


# is inside
def contains_properly_array(shape: ShapeType, arr: ShapeArray) -> NDArray[float]:
    shapely.prepare(shape.polygon)
    pass
