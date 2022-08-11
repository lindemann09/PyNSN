__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from .misc import numpy_array_2d
from typing import Union
from numpy.typing import ArrayLike, NDArray

# all functions are 2D arrays (at least) as fist arguments


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


def dist_coordinates(a_xy: ArrayLike,
                     b_xy: ArrayLike) -> NDArray:
    """Euclidean distances between coordinates (a and b)

    Note: At least one xy parameter has to a 2D array
    """
    d_xy = np.asarray(a_xy) - np.asarray(b_xy)  # type: ignore
    return np.hypot(d_xy[:, 0], d_xy[:, 1])


def dist_dots(a_xy: ArrayLike, a_diameter: ArrayLike,
              b_xy: ArrayLike, b_diameter: ArrayLike) -> NDArray:
    """Euclidean distances between dots (a and b)

    Note: At least one xy parameter has to a 2D array
    """
    # distance between centers minus the radii
    return dist_coordinates(a_xy, b_xy) \
        - (np.asarray(a_diameter) + np.asarray(b_diameter)) / 2


def xy_dist_rectangles(a_xy: ArrayLike, a_sizes: ArrayLike,
                       b_xy: ArrayLike, b_sizes: ArrayLike) -> NDArray:
    """return distances on both axes between rectangles.
    negative numbers indicate overlap of edges along that dimension.

    Note: At least one xy  parameter has to a 2D array
    """

    max_overlap_dist = (np.asarray(a_sizes) + np.asarray(b_sizes)) / 2
    diff = np.asarray(a_xy) - np.asarray(b_xy)  # type: ignore
    return np.abs(diff) - max_overlap_dist


def dist_rectangles(a_xy: ArrayLike, a_sizes: ArrayLike,
                    b_xy: ArrayLike, b_sizes: ArrayLike) -> NDArray:
    """Distance between rectangles

    negative distances indicate overlap and represent the
    size of the minimum overlap
    Note: At least one xy parameter has to a 2D array
    """
    d_xy = xy_dist_rectangles(a_xy, a_sizes, b_xy, b_sizes)
    # columns both negative -
    both_neg = np.flatnonzero(np.all(d_xy < 0, axis=1))
    # set make largest neg number positive for later dist calculation
    tmp = d_xy[both_neg, :]
    i = tmp[:, 0] < tmp[:, 1]
    tmp[i, 1] = np.abs(tmp[i, 1])  # second is positive
    tmp[~i, 0] = np.abs(tmp[~i, 0])  # first is positive
    d_xy[both_neg, :] = tmp
    # all other negatives are ignored and only one dimension is considered
    d_xy[np.where(d_xy < 0)] = 0
    # calc euclidean distance
    rtn = np.hypot(d_xy[:, 0], d_xy[:, 1])
    # set overlaps negative
    rtn[both_neg] = -1 * rtn[both_neg]
    return rtn


def center_of_positions(xy):
    """ calc center of all positions"""
    min_max = np.array((np.min(xy, axis=0), np.max(xy, axis=0)))
    return np.reshape(min_max[1, :] - np.diff(min_max, axis=0) / 2, 2)
