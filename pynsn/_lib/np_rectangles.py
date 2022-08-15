"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import NDArray

from .np_tools import CombinationMatrx
from .np_coordinates import cartesian2polar

# all functions are 2D arrays (at least) as fist arguments


def xy_distances(a_xy: NDArray, a_sizes: NDArray,
                 b_xy: NDArray, b_sizes: NDArray) -> NDArray:
    """return distances on both axes between rectangles edges.
    negative numbers indicate overlap of corners along that dimension.

    Note: At least one xy  parameter has to be a 2D array
    """

    return np.abs(a_xy - b_xy) - (a_sizes + b_sizes) / 2  # type: ignore


def overlaps(a_xy: NDArray, a_sizes: NDArray,
             b_xy: NDArray, b_sizes: NDArray,
             minimum_distance: float = 0) -> NDArray:
    """True if rectangles overlap

    Note: At least one xy parameter has to be a 2D array
    """
    # columns both negative -> it overlaps
    return np.all(xy_distances(a_xy, a_sizes, b_xy, b_sizes)
                  < minimum_distance, axis=1)


def overlap_matrix(xy: NDArray, sizes: NDArray,
                   minimum_distance: float = 0) -> NDArray:
    """Matrix with overlaps (True/False) between the rectangles"""
    mtx = CombinationMatrx(xy.shape[0])
    dist = overlaps(a_xy=xy[mtx.idx_a, :],
                    a_sizes=sizes[mtx.idx_a, :],
                    b_xy=xy[mtx.idx_b, :],
                    b_sizes=sizes[mtx.idx_b, :],
                    minimum_distance=minimum_distance)
    return mtx.fill(values=dist)


def distances(a_xy: NDArray, a_sizes: NDArray,
              b_xy: NDArray, b_sizes: NDArray) -> NDArray:
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


def distance_matrix(xy: NDArray, sizes: NDArray) -> NDArray:
    """Return matrix with distance between the rectangles"""
    xy = np.array(xy)
    size = np.array(sizes)
    mtx = CombinationMatrx(xy.shape[0])
    dist = distances(a_xy=xy[mtx.idx_a, :],
                     a_sizes=size[mtx.idx_a, :],
                     b_xy=xy[mtx.idx_b, :],
                     b_sizes=size[mtx.idx_b, :])
    return mtx.fill(values=dist)


def overlap_with_dots(rect_xy: NDArray, rect_sizes: NDArray,
                      dot_xy: NDArray, dot_diameter: NDArray) -> NDArray:
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
    return np.all(dist_lbrt < 0, axis=1)


def inside_rectangle(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates with rectangles are fully inside the rectangle/rectangles"""
    raise NotImplementedError()


def inside_dots(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates with rectangles are fully inside the dot/dots"""
    raise NotImplementedError()


def center_edge_distance(angle: NDArray, rect_sizes: NDArray) -> NDArray:
    """Distance between rectangle center and rectangle edge along the  line
    in direction of `angle`.
    """

    if rect_sizes.ndim == 1:
        rect_sizes = np.ones((angle.shape[0], 1)) * rect_sizes

    l_inside = np.empty(len(angle))
    # find vertical relations
    v_rel = (np.pi-np.pi/4 >= abs(angle)) & (abs(angle) > np.pi/4)
    # vertical relation: in case line cut rectangle at the top or bottom corner
    idx = np.flatnonzero(v_rel)
    l_inside[idx] = rect_sizes[idx, 1] / 2 * np.cos(np.pi/2 - angle[idx])
    # horizontal relation: in case line cut rectangle at the left or right corner
    idx = np.flatnonzero(~v_rel)
    l_inside[idx] = rect_sizes[idx, 0] / 2 * np.cos(angle[idx])

    return l_inside


def required_displacement(a_xy: NDArray, a_sizes: NDArray,
                          b_xy: NDArray, b_sizes: NDArray,
                          minimum_distance: float = 0) -> NDArray:
    """Minimum required displacement of object B to have no overlap with
    Object A

    If 0, there is no overlap.

    Movement direction is always along the line of the two object center
    """

    assert (a_xy.shape == a_sizes.shape ==
            b_xy.shape == b_sizes.shape)  # TODO could me more flecible

    xy_diff = a_xy - b_xy  # xy distance between object center
    min_xy_dist = ((a_sizes + b_sizes) / 2) + minimum_distance
    overlaps = np.all(np.abs(xy_diff) - min_xy_dist < 0, axis=1)

    rtn = np.zeros(xy_diff.shape)

    # only process overlapping items further
    i_over = np.flatnonzero(overlaps)
    xy_diff = xy_diff[i_over, :]
    min_xy_dist = min_xy_dist[i_over, :]

    # euclidean distances  and directions between center
    dist_direct = cartesian2polar(xy_diff[i_over, :])  # distance and direction

    # adapt distance to minimum required movement
    # dist(x)**2 + dist(y)**2 = euclidean_dist(xy)**2
    # set e.g. x to minimum dist:
    #   sqrt(min_dist(x)**2 + dist(y)**2) = min_euclidean_dist(xy)
    #
    # set the smaller xy distance to minimum and calc required movement along the
    # line between object center (euclidean distances)

    # find x smaller then y distance
    # if equals, randomly choose some and treat them like x_smaller
    x_smaller = np.abs(xy_diff[:, 0]) < np.abs(xy_diff[:, 1])

    # TODO
    # i_xy_equal = np.flatnonzero(xy_diff[:, 0] == xy_diff[:, 1])
    # if len(i_xy_equal) > 0:
    #    # randomly delete some items and add the remained to x_smaller
    #    r_del = np.random.choice([False, True], size=len(i_xy_equal))
    #    x_smaller = np.append(x_smaller, np.delete(i_xy_equal, r_del))

    # calc target distances
    target_eucl_dist = np.empty(len(i_over))
    # x smaller -> set x to min_dist
    i = np.flatnonzero(x_smaller)
    target_eucl_dist[i] = np.sqrt(min_xy_dist[i, 0] ** 2
                                  + xy_diff[i, 1] ** 2)
    # x larger -> set y to min_dist
    i = np.flatnonzero(~x_smaller)
    target_eucl_dist[i] = np.sqrt(xy_diff[i, 0] ** 2
                                  + min_xy_dist[i, 1] ** 2)

    # subtract distance inside rectangles
    l_inside_a = center_edge_distance(dist_direct[:, 1], a_sizes[i_over, :])
    l_inside_b = center_edge_distance(dist_direct[:, 1], b_sizes[i_over, :])
    dist_direct[:, 0] = target_eucl_dist - l_inside_a - l_inside_b

    # insert in return array
    rtn[i_over, :] = dist_direct

    return rtn
