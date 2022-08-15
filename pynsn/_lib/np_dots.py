"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import NDArray

from .np_tools import CombinationMatrx
from . import np_coordinates

# all functions are 2D arrays (at least) as fist arguments


def distances(a_xy: NDArray, a_diameter: NDArray,
              b_xy: NDArray, b_diameter: NDArray) -> NDArray:
    """Euclidean distances between dots (a and b)

    Note: At least one xy parameter has to be a 2D array
    """

    # distance between centers minus the radii
    object_expension = (np.asarray(a_diameter) +
                        np.asarray(b_diameter)) / 2
    d_xy = np.asarray(b_xy) - np.asarray(a_xy)  # type: ignore
    return np.hypot(d_xy[:, 0], d_xy[:, 1]) - object_expension


def inside_rectangles(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates with dots are fully inside the rectangle/rectangles"""
    raise NotImplementedError()


def inside_dots(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates with dots are fully inside the dot/dots"""
    raise NotImplementedError()


def distance_matrix(xy: NDArray, diameter: NDArray) -> NDArray:
    """Return matrix with distance between the dots"""
    xy = np.asarray(xy)
    diameter = np.asarray(diameter)
    mtx = CombinationMatrx(xy.shape[0])

    object_expension = (diameter[mtx.idx_b] + diameter[mtx.idx_a]) / 2
    d_xy = xy[mtx.idx_b, :] - xy[mtx.idx_a, :]  # type: ignore
    dist = np.hypot(d_xy[:, 0], d_xy[:, 1]) - object_expension

    return mtx.fill(values=dist)
