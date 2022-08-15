"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import NDArray

from .np_tools import CombinationMatrx, as_vector, make_array2d
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


def center_edge_distance(angle: NDArray, rect_sizes: NDArray) -> NDArray[np.floating]:
    """Distance between rectangle center and rectangle edge along the  line
    in direction of `angle`.
    """

    if rect_sizes.ndim == 1:
        rect_sizes = np.ones((angle.shape[0], 1)) * rect_sizes

    l_inside = np.empty(len(angle))
    # find vertical relations
    v_rel = (np.pi-np.pi/4 >= abs(angle)) & (abs(angle) > np.pi/4)
    # vertical relation: in case line cut rectangle at the top or bottom corner
    i = np.flatnonzero(v_rel)
    l_inside[i] = rect_sizes[i, 1] / 2 * np.cos(np.pi/2 - angle[i])
    # horizontal relation: in case line cut rectangle at the left or right corner
    i = np.flatnonzero(~v_rel)
    l_inside[i] = rect_sizes[i, 0] / 2 * np.cos(angle[i])

    return l_inside


GO ON HERE: TEST THESE FUNCTIONS


def required_displacement_with_rects(rects_a_xy: NDArray,
                                     rects_a_sizes: NDArray,
                                     rects_b_xy: NDArray,
                                     rects_b_sizes: NDArray,
                                     minimum_distance: float = 0) -> NDArray:
    """Minimum required displacement of rectangle B to have no overlap with
     rectangle A. Array with replacement vectors in polar coordinates.

    Two objects do not overlap if replacements[:, 0] == 0.

    Returns:
        replacements as array of vectors in polar coordinates
        if distance (replacements[:, 0]) is 0, no replacement required.

        Calculated movement direction (replacements[:, 0]) is always along the
        line between the two object center.
    """

    #
    # all arrays are 2D and have the same n_rows
    #
    # fjnd two 2d array with more than one row
    n_rows = 1
    for arr in (rects_a_xy, rects_a_sizes, rects_b_xy, rects_b_sizes):
        if arr.ndim > 1 and arr.shape[0] > 1:
            n_rows = arr.shape[0]
            break
    rects_a_xy = make_array2d(rects_a_xy, n_rows=n_rows)
    rects_a_sizes = make_array2d(rects_a_sizes, n_rows=n_rows)
    rects_b_xy = make_array2d(rects_b_xy, n_rows=n_rows)
    rects_b_sizes = make_array2d(rects_b_sizes, n_rows=n_rows)

    #
    # find overlapping objects
    #
    xy_dist = rects_b_xy - rects_a_xy  # type: ignore
    min_xy_dist = ((rects_a_sizes + rects_b_sizes) / 2) + minimum_distance
    overlapping = np.all(np.abs(xy_dist) - min_xy_dist < 0, axis=1)
    i_over = np.flatnonzero(overlapping)

    #
    # calc dist and direction for overlapping objects
    #
    rtn = np.zeros(xy_dist.shape)

    # only process overlapping items further
    min_spatial_rel = _min_spatial_relation(xy_dist=xy_dist[i_over, :],
                                            min_xy_dist=min_xy_dist[i_over, :])
    # subtract distance inside rectangles from euclidean distance
    l_inside_a = center_edge_distance(min_spatial_rel[:, 1],
                                      rects_a_sizes[i_over, :])
    l_inside_b = center_edge_distance(min_spatial_rel[:, 1],
                                      rects_b_sizes[i_over, :])
    min_spatial_rel[:, 0] = min_spatial_rel[:, 0] - (l_inside_a + l_inside_b)

    # insert calculated replacements into return array
    rtn[i_over, :] = min_spatial_rel

    return rtn


def required_displacement_with_dot(rects_xy: NDArray,
                                   rects_sizes: NDArray,
                                   dots_xy: NDArray,
                                   dots_diameter: NDArray,
                                   minimum_distance: float = 0) -> NDArray:
    """Minimum required displacement of rectangle B to have no overlap with
     rectangle A. Array with replacement vectors in polar coordinates.

    Two objects do not overlap if replacements[:, 0] == 0.

    Returns:
        replacements as array of vectors in polar coordinates
        if distance (replacements[:, 0]) is 0, no replacement required.

        Calculated movement direction (replacements[:, 0]) is always along the
        line between the two object center.
    """

    #
    # all arrays are 2D and have the same n_rows
    #
    # fjnd two 2d array with more than one row
    n_rows = 1
    for arr in (rects_xy, rects_sizes, dots_xy):
        if arr.ndim > 1 and arr.shape[0] > 1:
            n_rows = arr.shape[0]
            break
    rects_xy = make_array2d(rects_xy, n_rows=n_rows)
    rects_sizes = make_array2d(rects_sizes, n_rows=n_rows)
    dots_xy = make_array2d(dots_xy, n_rows=n_rows)
    dots_diameter = as_vector(dots_diameter)
    dots_rect_sizes = np.transpose(
        np.ones((2, dots_diameter.shape[0])) * dots_diameter)

    #
    # find xy overlapping objects (simplified: treat dot as rect)
    #
    xy_dist = dots_xy - rects_xy  # type: ignore
    min_xy_dist = ((rects_sizes + dots_rect_sizes) / 2) + minimum_distance
    overlapping = np.all(np.abs(xy_dist) - min_xy_dist < 0, axis=1)
    i_over = np.flatnonzero(overlapping)

    #
    # calc dist and direction for overlapping objects
    #
    rtn = np.zeros(xy_dist.shape)

    # only process overlapping items further
    min_spatial_rel = _min_spatial_relation(xy_dist=xy_dist[i_over, :],
                                            min_xy_dist=min_xy_dist[i_over, :])
    # subtract distance inside rectangles from euclidean distance
    l_inside_rect = center_edge_distance(min_spatial_rel[:, 1],
                                         rects_sizes[i_over, :])
    l_inside_dot = dots_diameter/2
    min_spatial_rel[:, 0] = min_spatial_rel[:, 0] \
        - (l_inside_rect + l_inside_dot)
    # FIXME check if replacement is at all required here. spatial relation might be smaller than inside_rect + inside dot
    # insert calculated replacements into return array
    rtn[i_over, :] = min_spatial_rel

    return rtn


def _min_spatial_relation(xy_dist: NDArray, min_xy_dist: NDArray) -> NDArray[np.floating]:
    """Calculate required minimum Euclidean distance and direction
    (polar coordinates) from xy distances and min_minimum required xy distances

    dist(x)**2 + dist(y)**2 = euclidean_dist(xy)**2
    set e.g. x to minimum dist:
      sqrt(min_dist(x)**2 + dist(y)**2) = min_euclidean_dist(xy)

    Set the smaller xy distance to minimum dist and calc required movement
    long the line between object centers (Euclidean distances minus length of
    line inside objects).
    """

    # find x smaller then y distance
    # if equals, randomly choose some and treat them like x_smaller
    x_smaller = np.abs(xy_dist[:, 0]) < np.abs(xy_dist[:, 1])
    xy_equal = xy_dist[:, 0] == xy_dist[:, 1]
    if np.sum(xy_equal) > 0:
        # randomly take some (via filtering) and add to x_smaller
        filter_vector = np.random.choice([False, True], size=len(xy_equal))
        x_smaller = x_smaller | (xy_equal & filter_vector)

    # calc required euclidean distances between object center
    min_eucl_dist = np.empty(xy_dist.shape[0])
    i = np.flatnonzero(x_smaller)  # x smaller -> set x to min_xy_dist
    min_eucl_dist[i] = np.hypot(min_xy_dist[i, 0], xy_dist[i, 1])
    i = np.flatnonzero(~x_smaller)  # x larger -> set y to min_xy_dist
    min_eucl_dist[i] = np.hypot(xy_dist[i, 0], min_xy_dist[i, 1])

    # return array of polar vectors
    # [[euclidean distances, directions between center]]
    directions = np.arctan2(xy_dist[:, 1], xy_dist[:, 0])
    return np.array([min_eucl_dist, directions]).T
