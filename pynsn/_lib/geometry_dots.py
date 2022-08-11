"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import geometry_coordinates as geo_coord
from .misc import CombinationMatrx

# all functions are 2D arrays (at least) as fist arguments


@staticmethod
def distances(a_xy: ArrayLike, a_diameter: ArrayLike,
              b_xy: ArrayLike, b_diameter: ArrayLike) -> NDArray:
    """Euclidean distances between dots (a and b)

    Note: At least one xy parameter has to be a 2D array
    """
    # distance between centers minus the radii
    object_expension = (np.asarray(a_diameter) +
                        np.asarray(b_diameter)) / 2
    return geo_coord.distances(a_xy, b_xy) - object_expension


@staticmethod
def inside_rectangles(xy: ArrayLike, sizes: ArrayLike) -> NDArray:
    """bool array indicates with dots are fully inside the rectangle/rectangles"""
    raise NotImplementedError()


@staticmethod
def inside_dots(xy: ArrayLike, sizes: ArrayLike) -> NDArray:
    """bool array indicates with dots are fully inside the dot/dots"""
    raise NotImplementedError()


@staticmethod
def distance_matrix(xy: ArrayLike, diameter: ArrayLike) -> NDArray:
    """Return matrix with distance between the dots"""
    xy = np.array(xy)
    diameter = np.array(diameter)
    mtx = CombinationMatrx(xy.shape[0])
    dist = distances(a_xy=xy[mtx.idx_a, :],
                     a_diameter=diameter[mtx.idx_a],
                     b_xy=xy[mtx.idx_b, :],
                     b_diameter=diameter[mtx.idx_b])
    return mtx.fill(values=dist)
