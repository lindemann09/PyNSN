__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
import itertools
from typing import Union
from numpy.typing import ArrayLike, NDArray

# all functions are 2D arrays (at least) as fist arguments


class _CombinationMatrx(object):
    """Symmetric combination matrix
    helper class"""

    def __init__(self, n_items: int) -> None:
        idx = np.array(list(itertools.combinations(range(n_items), r=2)))
        self.idx_a = idx[:, 0]
        self.idx_b = idx[:, 1]
        self.matrix = np.full((n_items, n_items), np.nan)

    @property
    def n_combinations(self) -> float:
        return len(self.idx_a)

    def fill(self, values) -> NDArray:
        """returns combination matrix with values

        """
        self.matrix[self.idx_a, self.idx_b] = values
        self.matrix[self.idx_b, self.idx_a] = values
        return self.matrix


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


def coordinates_distances(a_xy: ArrayLike,
                          b_xy: ArrayLike) -> NDArray:
    """Euclidean distances between coordinates (a and b)

    Note: At least one xy parameter has to be a 2D array
    """
    d_xy = np.asarray(a_xy) - np.asarray(b_xy)  # type: ignore
    return np.hypot(d_xy[:, 0], d_xy[:, 1])


def dots_distances(a_xy: ArrayLike, a_diameter: ArrayLike,
                   b_xy: ArrayLike, b_diameter: ArrayLike) -> NDArray:
    """Euclidean distances between dots (a and b)

    Note: At least one xy parameter has to be a 2D array
    """
    # distance between centers minus the radii
    object_expension = (np.asarray(a_diameter) + np.asarray(b_diameter)) / 2
    return coordinates_distances(a_xy, b_xy) - object_expension


def rectangles_xy_distances(a_xy: ArrayLike, a_sizes: ArrayLike,
                            b_xy: ArrayLike, b_sizes: ArrayLike) -> NDArray:
    """return distances on both axes between rectangles.
    negative numbers indicate overlap of edges along that dimension.

    Note: At least one xy  parameter has to be a 2D array
    """

    object_expension = (np.asarray(a_sizes) + np.asarray(b_sizes)) / 2
    diff = np.asarray(a_xy) - np.asarray(b_xy)  # type: ignore
    return np.abs(diff) - object_expension


def rectangles_distances(a_xy: ArrayLike, a_sizes: ArrayLike,
                         b_xy: ArrayLike, b_sizes: ArrayLike) -> NDArray:
    """Distance between rectangles

    negative distances indicate overlap and represent the
    size of the minimum overlap
    Note: At least one xy parameter has to be a 2D array
    """
    d_xy = rectangles_xy_distances(a_xy, a_sizes, b_xy, b_sizes)
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


def rectangles_overlap(a_xy: ArrayLike, a_sizes: ArrayLike,
                       b_xy: ArrayLike, b_sizes: ArrayLike) -> NDArray:
    """True if rectangles overlap

    Note: At least one xy parameter has to be a 2D array
    """
    d_xy = rectangles_xy_distances(a_xy, a_sizes, b_xy, b_sizes)
    # columns both negative -
    return np.all(d_xy < 0, axis=1)


def rectangles_distance_matrix(xy: ArrayLike, sizes: ArrayLike) -> NDArray:
    """Return matrix with distance between the rectangles"""
    xy = np.array(xy)
    size = np.array(sizes)
    mtx = _CombinationMatrx(xy.shape[0])
    dist = rectangles_distances(a_xy=xy[mtx.idx_a, :],
                                a_sizes=size[mtx.idx_a, :],
                                b_xy=xy[mtx.idx_b, :],
                                b_sizes=size[mtx.idx_b, :])
    return mtx.fill(values=dist)


def rectangles_overlap_matrix(xy: ArrayLike, sizes: ArrayLike) -> NDArray:
    """Return matrix with overlaps (True/False) between the rectangles"""
    xy = np.array(xy)
    size = np.array(sizes)
    mtx = _CombinationMatrx(xy.shape[0])
    dist = rectangles_overlap(a_xy=xy[mtx.idx_a, :],
                              a_sizes=size[mtx.idx_a, :],
                              b_xy=xy[mtx.idx_b, :],
                              b_sizes=size[mtx.idx_b, :])
    return mtx.fill(values=dist)


def dots_distance_matrix(xy: ArrayLike, diameter: ArrayLike) -> NDArray:
    """Return matrix with distance between the dots"""
    xy = np.array(xy)
    diameter = np.array(diameter)
    mtx = _CombinationMatrx(xy.shape[0])
    dist = dots_distances(a_xy=xy[mtx.idx_a, :],
                          a_diameter=diameter[mtx.idx_a],
                          b_xy=xy[mtx.idx_b, :],
                          b_diameter=diameter[mtx.idx_b])
    return mtx.fill(values=dist)


def coordinates_distance_matrix(xy: ArrayLike) -> NDArray:
    """Return matrix with distance between the dots"""
    xy = np.array(xy)
    n = xy.shape[0]
    mtx = _CombinationMatrx(xy.shape[0])
    dist = coordinates_distances(a_xy=xy[mtx.idx_a, :],
                                 b_xy=xy[mtx.idx_b, :])
    return mtx.fill(values=dist)


def center_of_positions(xy):
    """ calc center of all positions"""
    min_max = np.array((np.min(xy, axis=0), np.max(xy, axis=0)))
    return np.reshape(min_max[1, :] - np.diff(min_max, axis=0) / 2, 2)
