"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .misc import CombinationMatrx

# all functions are 2D arrays (at least) as fist arguments


def distance_matrix(xy: ArrayLike) -> NDArray:
    """Return matrix with distance between the dots"""
    xy = np.array(xy)
    n = xy.shape[0]
    mtx = CombinationMatrx(xy.shape[0])
    dist = distances(a_xy=xy[mtx.idx_a, :],
                     b_xy=xy[mtx.idx_b, :])
    return mtx.fill(values=dist)


def center_of_coordinates(xy: ArrayLike) -> NDArray:
    """ calc center of all positions"""
    min_max = np.array((np.min(xy, axis=0), np.max(xy, axis=0)))
    return np.reshape(min_max[1, :] - np.diff(min_max, axis=0) / 2, 2)


def polar2cartesian(polar: ArrayLike) -> NDArray:
    """polar has to be an 2D-array representing polar coordinates (radius, angle)"""
    polar = np.asarray(polar)
    return np.array([polar[:, 0] * np.cos(polar[:, 1]),
                    polar[:, 0] * np.sin(polar[:, 1])]).T


def cartesian2polar(xy: ArrayLike,
                    radii_only: bool = False) -> NDArray:
    """polar coordinates (radius, angle)

    if only radii required you may consider radii_only=True for faster
    processing

    xy has to be a 2D array
    """
    xy = np.asarray(xy)
    rtn = np.hypot(xy[:, 0], xy[:, 1])
    if not radii_only:
        # add angle column
        rtn = np.array([rtn, np.arctan2(xy[:, 1], xy[:, 0])]).T
    return rtn


def cartesian2image_coordinates(xy: ArrayLike,
                                image_size: ArrayLike) -> NDArray:
    """convert cartesian to image coordinates with (0,0) at top left and
    reversed y axis

    xy has to be a 2D array

    """
    return (np.asarray(xy) * [1, -1]) + np.asarray(image_size) / 2


def distances(a_xy: ArrayLike,
              b_xy: ArrayLike) -> NDArray:
    """Euclidean distances between coordinates (a and b)

    Note: At least one xy parameter has to be a 2D array
    """
    d_xy = np.asarray(a_xy) - np.asarray(b_xy)  # type: ignore
    return np.hypot(d_xy[:, 0], d_xy[:, 1])
