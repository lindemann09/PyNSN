from typing import Callable
import numpy as np
import shapely
from .shape_array import ShapeArray


def _spatial_relations(shape_array: ShapeArray,
                       relations_fnc: Callable):
    """helper function returning the relation between polygons
    """
    l = len(shape_array.polygons)
    rtn = np.full((l, l), np.nan)
    idx = list(range(l))
    while len(idx):
        r = idx.pop()  # starts with last element (i.e. last row)
        rtn[r, idx] = relations_fnc(shape_array.polygons[idx],
                                    shape_array.polygons[r])
    return rtn


def dwithin(shape_array: ShapeArray, dist: float):
    return _spatial_relations(shape_array,
                              lambda a, b: shapely.dwithin(a, b, distance=dist))


def distance(shape_array: ShapeArray):
    return _spatial_relations(shape_array, shapely.distance)


def intersects(shape_array: ShapeArray):
    return _spatial_relations(shape_array, shapely.intersects)


def contains_properly(shape_array: ShapeArray):
    return _spatial_relations(shape_array, shapely.contains_properly)
