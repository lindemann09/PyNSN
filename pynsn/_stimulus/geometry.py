from typing import Callable

import numpy as _np
import shapely as _shapely
from numpy.typing import NDArray

from .. import _stimulus
from .._shapes.circ_shapes import (distance_circ_circ, ellipse_diameter,
                                  ellipse_perimeter, is_circ_in_circ)
from .._shapes.shapes import is_in_shape
from .shape_array import (ShapeArray, distance_circ_dot_array,
                                    distance_circ_ellipse_array)


def _matrix_spatial_relations(shape_array: ShapeArray,
                              relations_fnc: Callable):
    """helper function returning the relation between polygons
    """
    l = len(shape_array.polygons)
    rtn = _np.full((l, l), _np.nan)
    idx = list(range(l))
    while len(idx):
        r = idx.pop()  # starts with last element (i.e. last row)
        rtn[r, idx] = relations_fnc(shape_array.polygons[idx],
                                    shape_array.polygons[r])
    return rtn


def matrix_dwithin(shape_array: _stimulus.ShapeArray, dist: float):
    return _matrix_spatial_relations(shape_array,
                                     lambda a, b: _shapely.dwithin(a, b, distance=dist))


def matrix_distance(shape_array: _stimulus.ShapeArray):
    return _matrix_spatial_relations(shape_array, _shapely.distance)


def matrix_intersects(shape_array: _stimulus.ShapeArray):
    return _matrix_spatial_relations(shape_array, _shapely.intersects)


def matrix_isinside(shape_array: _stimulus.ShapeArray):
    return _matrix_spatial_relations(shape_array, _shapely.contains_properly)
