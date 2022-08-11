"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import geometry_coordinates as geo_coord
from .misc import CombinationMatrx

# all functions are 2D arrays (at least) as fist arguments


def xy_distances(a_xy: ArrayLike, a_sizes: ArrayLike,
                 b_xy: ArrayLike, b_sizes: ArrayLike) -> NDArray:
    """return distances on both axes between rectangles.
    negative numbers indicate overlap of edges along that dimension.

    Note: At least one xy  parameter has to be a 2D array
    """

    object_expension = (np.asarray(a_sizes) + np.asarray(b_sizes)) / 2
    diff = np.asarray(a_xy) - np.asarray(b_xy)  # type: ignore
    return np.abs(diff) - object_expension


def distances(a_xy: ArrayLike, a_sizes: ArrayLike,
              b_xy: ArrayLike, b_sizes: ArrayLike) -> NDArray:
    """Distance between rectangles

    negative distances indicate overlap and represent the
    size of the minimum overlap
    Note: At least one xy parameter has to be a 2D array
    """
    d_xy = xy_distances(a_xy, a_sizes, b_xy, b_sizes)
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


def overlaps(a_xy: ArrayLike, a_sizes: ArrayLike,
             b_xy: ArrayLike, b_sizes: ArrayLike) -> NDArray:
    """True if rectangles overlap

    Note: At least one xy parameter has to be a 2D array
    """
    d_xy = xy_distances(a_xy, a_sizes, b_xy, b_sizes)
    # columns both negative -
    return np.all(d_xy < 0, axis=1)


def distance_matrix(xy: ArrayLike, sizes: ArrayLike) -> NDArray:
    """Return matrix with distance between the rectangles"""
    xy = np.array(xy)
    size = np.array(sizes)
    mtx = CombinationMatrx(xy.shape[0])
    dist = distances(a_xy=xy[mtx.idx_a, :],
                     a_sizes=size[mtx.idx_a, :],
                     b_xy=xy[mtx.idx_b, :],
                     b_sizes=size[mtx.idx_b, :])
    return mtx.fill(values=dist)


def overlap_matrix(xy: ArrayLike, sizes: ArrayLike) -> NDArray:
    """Matrix with overlaps (True/False) between the rectangles"""
    xy = np.array(xy)
    size = np.array(sizes)
    mtx = CombinationMatrx(xy.shape[0])
    dist = overlaps(a_xy=xy[mtx.idx_a, :],
                    a_sizes=size[mtx.idx_a, :],
                    b_xy=xy[mtx.idx_b, :],
                    b_sizes=size[mtx.idx_b, :])
    return mtx.fill(values=dist)


def overlap_with_dots(rect_xy: ArrayLike, rect_sizes: ArrayLike,
                      dot_xy: ArrayLike, dot_diameter: ArrayLike) -> NDArray:
    """Array with overlaps (True, False) of the rectangle and the dot/dots"""
    rect_xy = np.asarray(rect_xy)
    rect_sizes2 = np.asarray(rect_sizes) / 2
    dot_xy = np.asarray(dot_xy)
    dot_diameter = np.asarray(dot_diameter)
    dot_diameter = dot_diameter.reshape((dot_diameter.shape[0], 1))

    # positive dist -> dot outside rect, negative is a potential overlap
    dist_left_bottom = (rect_xy - rect_sizes2) - dot_xy \
        - dot_diameter  # +y = dot below,  +x = dot left
    dist_right_top = dot_xy - (rect_xy + rect_sizes2) \
        - dot_diameter  # +y = dot top, +x dot right

    dist_lbrt = np.append(dist_left_bottom, dist_right_top, axis=1)
    # if dot is overlapping -> all distances negative
    overlaps = np.all(dist_lbrt < 0, axis=1)
    return overlaps


def inside_rectangle(xy: ArrayLike, sizes: ArrayLike) -> NDArray:
    """bool array indicates with rectangles are fully inside the rectangle/rectangles"""
    raise NotImplementedError()


def inside_dots(xy: ArrayLike, sizes: ArrayLike) -> NDArray:
    """bool array indicates with rectangles are fully inside the dot/dots"""
    raise NotImplementedError()


WIETER CHECK ALL gemortry imports
