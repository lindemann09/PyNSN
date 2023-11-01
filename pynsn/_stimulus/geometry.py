from typing import Callable, Union

import numpy as np
import shapely
from numpy.typing import NDArray

from .. import _shapes
from .. import _stimulus
from .._shapes.geometry import ellipse_diameter


def distance_circ_dot_array(obj: Union[_shapes.Point2D, _shapes.CircularShapeType],
                            dots_xy: NDArray[np.float_],
                            dots_diameter: NDArray[np.float_]) -> NDArray[np.float_]:
    """Distances circular shape or Point to multiple dots
    """
    d_xy = dots_xy - obj.xy
    if isinstance(obj, _shapes.Point2D):
        circ_dia = 0
    elif isinstance(obj, _shapes.Dot):
        circ_dia = obj.diameter
    elif isinstance(obj, _shapes.Ellipse):
        circ_dia = ellipse_diameter(
            size=np.atleast_2d(obj.size),
            theta=np.arctan2(d_xy[:, 1], d_xy[:, 0]))  # ellipse radius to each dot in the array
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(obj)}")

    # center dist - radius_a - radius_b
    return np.hypot(d_xy[:, 0], d_xy[:, 1]) - (dots_diameter + circ_dia) / 2


def distance_circ_ellipse_array(obj: Union[_shapes.Point2D, _shapes.CircularShapeType],
                                ellipses_xy: NDArray[np.float_],
                                ellipse_sizes: NDArray[np.float_]) -> NDArray[np.float_]:
    """Distance circular shape or Point2D to multiple ellipses
    """
    d_xy = ellipses_xy - obj.xy
    theta = np.arctan2(d_xy[:, 0], d_xy[:, 1])
    # radii of ellipses in array to circ_shape
    ellipse_dia = ellipse_diameter(size=ellipse_sizes, theta=theta)
    if isinstance(obj, _shapes.Point2D):
        shape_dia = 0
    elif isinstance(obj, _shapes.Dot):
        shape_dia = obj.diameter
    elif isinstance(obj, _shapes.Ellipse):
        shape_dia = ellipse_diameter(size=np.atleast_2d(obj.size),
                                     theta=theta)  # ellipse radius to each ellipse in the array
    else:
        raise RuntimeError(f"Unknown circular shape type: {type(obj)}")

    # center dist - radius_a - radius_b
    return np.hypot(d_xy[:, 0], d_xy[:, 1]) - (ellipse_dia + shape_dia)/2


def _matrix_spatial_relations(shape_array: _stimulus.ShapeArray,
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


def matrix_dwithin(shape_array: _stimulus.ShapeArray, dist: float):
    return _matrix_spatial_relations(shape_array,
                                     lambda a, b: shapely.dwithin(a, b, distance=dist))


def matrix_distance(shape_array: _stimulus.ShapeArray):
    return _matrix_spatial_relations(shape_array, shapely.distance)


def matrix_intersects(shape_array: _stimulus.ShapeArray):
    return _matrix_spatial_relations(shape_array, shapely.intersects)


def matrix_isinside(shape_array: _stimulus.ShapeArray):
    return _matrix_spatial_relations(shape_array, shapely.contains_properly)
