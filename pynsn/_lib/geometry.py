"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import NDArray, ArrayLike

from . import rng

# all functions are 2D arrays (at least) as fist arguments


def center_of_coordinates(xy: NDArray) -> NDArray:
    """ calc center of all positions"""
    min_max = np.array((np.min(xy, axis=0), np.max(xy, axis=0)))
    return np.reshape(min_max[1, :] - np.diff(min_max, axis=0) / 2, 2)


def polar2cartesian(polar: ArrayLike) -> NDArray:
    """polar has to be an 2D-array representing polar coordinates (radius, angle)"""
    polar = np.atleast_2d(polar)
    return np.array([polar[:, 0] * np.cos(polar[:, 1]),
                    polar[:, 0] * np.sin(polar[:, 1])]).T


def cartesian2polar(xy: ArrayLike,
                    radii_only: bool = False) -> NDArray[np.floating]:
    """polar coordinates (radius, angle)

    if only radii required you may consider radii_only=True for faster
    processing

    xy has to be a 2D array
    """
    xy = np.atleast_2d(xy)
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


def jitter_identical_coordinates(xy: NDArray, jitter_size: float = 0.1) -> NDArray:
    """jitters points with identical position"""

    for idx, ref_object in enumerate(xy):
        # find identical positions
        identical = np.flatnonzero(np.all(np.equal(xy, ref_object), axis=1))
        if len(identical) > 1:
            for x in identical:  # jitter all identical positions
                if x != idx:
                    tmp = np.atleast_2d((jitter_size,
                                        rng.generator.random() * 2 * np.pi))
                    xy[x, :] = xy[x, :] - polar2cartesian(tmp)[0]

    return xy


def corner_tensor(rect_xy: NDArray, rect_sizes: NDArray) -> NDArray:
    """returns a tensor (n, 2, 4) with xy values of the four corners"""
    assert rect_xy.shape == rect_sizes.shape
    rtn = np.empty((rect_xy.shape[0], 2, 4))

    rect_sizes2 = rect_sizes / 2
    rtn[:, :, 1] = rect_xy + rect_sizes2  # right top
    rtn[:, :, 3] = rect_xy - rect_sizes2  # left bottom
    # left, top
    rtn[:, 0, 0] = rtn[:, 0, 3]
    rtn[:, 1, 0] = rtn[:, 1, 1]
    # right, bottom
    rtn[:, 0, 2] = rtn[:, 0, 1]
    rtn[:, 1, 2] = rtn[:, 1, 3]
    return rtn


def dots_outer_points(dot_xy: NDArray, dot_radii: NDArray) -> NDArray:
    """returns a tensor (n, 2, 4) with four outer points (left, top, right, bottom)
    on the cardinal axes"""
    assert dot_xy.shape[0] == dot_radii.shape[0]

    # make tensor (n, 2, 4) with dot_xy at each [:, :, x]
    rtn = dot_xy.reshape((dot_xy.shape[0], dot_xy.shape[1], 1)) * \
        np.ones((dot_xy.shape[0], 2, 4))
    rtn[:, 0, 0] = rtn[:, 0, 0] - dot_radii  # left
    rtn[:, 1, 1] = rtn[:, 1, 1] + dot_radii  # top
    rtn[:, 0, 2] = rtn[:, 0, 2] + dot_radii  # right
    rtn[:, 1, 3] = rtn[:, 1, 3] - dot_radii  # bottom
    return rtn
