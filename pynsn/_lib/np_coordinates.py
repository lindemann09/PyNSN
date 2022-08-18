"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import NDArray, ArrayLike
from scipy.spatial.distance import cdist

from .._lib import rng

from .np_tools import CombinationMatrx, as_array2d

# all functions are 2D arrays (at least) as fist arguments


def distance_matrix(xy: NDArray) -> NDArray:
    """Return matrix with distance between the coordinates"""
    mtx = CombinationMatrx(xy.shape[0])
    d_xy = xy[mtx.idx_b, :] - xy[mtx.idx_a, :]  # type: ignore
    dist = np.hypot(d_xy[:, 0], d_xy[:, 1])
    return mtx.fill(values=dist)


def center_of_coordinates(xy: NDArray) -> NDArray:
    """ calc center of all positions"""
    min_max = np.array((np.min(xy, axis=0), np.max(xy, axis=0)))
    return np.reshape(min_max[1, :] - np.diff(min_max, axis=0) / 2, 2)


def polar2cartesian(polar: NDArray) -> NDArray:
    """polar has to be an 2D-array representing polar coordinates (radius, angle)"""
    polar = np.asarray(polar)
    return np.array([polar[:, 0] * np.cos(polar[:, 1]),
                    polar[:, 0] * np.sin(polar[:, 1])]).T


def cartesian2polar(xy: ArrayLike,
                    radii_only: bool = False) -> NDArray[np.floating]:
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


def spatial_relation(a_xy, b_xy):
    """distance and angle between two coordinates"""
    # distance
    return cartesian2polar(b_xy - a_xy)  # type: ignore


def distances(a_xy: NDArray, b_xy: NDArray) -> NDArray:
    """Euclidean distances between coordinates (a and b)

    Note: At least one xy parameter has to be a 2D array
    """
    d_xy = b_xy - a_xy  # type: ignore
    return np.hypot(d_xy[:, 0], d_xy[:, 1])


def jitter_identical_coordinates(xy: NDArray, jitter_size: float = 0.1) -> NDArray:
    """jitters points with identical position"""

    for idx, ref_object in enumerate(xy):
        # find identical positions
        identical = np.flatnonzero(np.all(np.equal(xy, ref_object), axis=1))
        if len(identical) > 1:
            for x in identical:  # jitter all identical positions
                if x != idx:
                    r = as_array2d((jitter_size,
                                    rng.generator.random() * 2 * np.pi))
                    xy[x, :] = xy[x, :] - polar2cartesian(r)[0]

    return xy
