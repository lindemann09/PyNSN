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


def outer_points(dot_xy: NDArray, dot_diameter: NDArray) -> NDArray:
    """returns a tensor (n, 2, 4) with four outer points (left, top, right, bottom)
    on the cardinal axes"""
    assert dot_xy.shape[0] == dot_diameter.shape[0]

    # radius vector
    radii = dot_diameter.reshape((dot_diameter.shape[0], )) / 2

    # make tensor (n, 2, 4) with dot_xy at each [:, :, x]
    rtn = dot_xy.reshape((dot_xy.shape[0], dot_xy.shape[1], 1)) * \
        np.ones((dot_xy.shape[0], 2, 4))
    rtn[:, 0, 0] = rtn[:, 0, 0] - radii  # left
    rtn[:, 1, 1] = rtn[:, 1, 1] + radii  # top
    rtn[:, 0, 2] = rtn[:, 0, 2] + radii  # right
    rtn[:, 1, 3] = rtn[:, 1, 3] - radii  # bottom
    return rtn
