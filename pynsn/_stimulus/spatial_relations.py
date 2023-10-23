"""
calculating spatial relations optimized for Dots

that is, for spatial relation between dots shapely is not used and relations are
calculated based on diameter and position
"""

from numpy.typing import NDArray
import numpy as np
import shapely
from .shape_array import ShapeArray, ShapeType
from .shapes import Dot, Rectangle, Picture


def distance(a:ShapeType, b:ShapeType) -> float:
    if isinstance(a, Dot) and isinstance(b, Dot):
        # dot - dot
        return np.hypot(b.xy[0] - a.xy[0], b.xy[1] - a.xy[1]) \
                - (a.diameter + b.diameter)/2 # -> - radius_a - radius_b
    else:
        return shapely.distance(a.polygon, b.polygon)


def dwithin(a:ShapeType, b:ShapeType, dist:float) -> bool:
    if isinstance(a, Dot) and isinstance(b, Dot):
        # dot - dot
        return distance(a,b) <= dist
    else:
        return shapely.dwithin(a.polygon, b.polygon, distance=dist)

def contains_properly(a:ShapeType, b:ShapeType) -> bool: # is inside
    if isinstance(a, Dot) and isinstance(b, Dot):
        # dot - dot
        return distance(a, b) < -2 * b.diameter
    else:
        return shapely.contains_properly(a.polygon, b.polygon)

def distance_array(shape:ShapeType, arr:ShapeArray) -> NDArray[float]:
    shapely.prepare(shape.polygon)

    dia = arr.dot_diameter[arr.dot_ids]
    xy = arr.xy[arr.dot_ids, :]


def dwithin_array(shape:ShapeType, arr:ShapeArray, dist:float) -> NDArray[bool]:
    shapely.prepare(shape.polygon)
    pass

def contains_properly_array(shape:ShapeType, arr:ShapeArray) -> NDArray[float]: # is inside
    shapely.prepare(shape.polygon)
    pass