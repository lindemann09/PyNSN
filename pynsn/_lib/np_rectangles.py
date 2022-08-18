"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from types import NotImplementedType
from typing import Optional, Union
import numpy as np
from numpy.typing import NDArray

from .np_tools import CombinationMatrx, as_vector, as_array2d_nrow
from .np_coordinates import cartesian2polar
from . import np_dots
# all functions are 2D arrays (at least) as fist arguments


def overlap(a_xy: NDArray, a_sizes: NDArray,
            b_xy: NDArray, b_sizes: NDArray,
            minimum_distance: float = 0) -> Union[NDArray[np.bool_], np.bool_]:
    """True if rectangles overlap

    Note: At least one xy parameter has to be a 2D array
    """
    # columns both negative -> it overlaps
    xy_dist = np.abs(a_xy - b_xy) - (a_sizes + b_sizes) / 2  # type: ignore
    comp = xy_dist < minimum_distance
    if comp.shape[0] > 1:
        return np.all(comp, axis=1)
    else:
        return np.all(comp)


def overlap_matrix(xy: NDArray, sizes: NDArray,
                   minimum_distance: float = 0) -> NDArray:
    """Matrix with overlaps (True/False) between the rectangles"""
    assert xy.shape == sizes.shape
    mtx = CombinationMatrx(xy.shape[0])
    dist = overlap(a_xy=xy[mtx.idx_a, :],
                   a_sizes=sizes[mtx.idx_a, :],
                   b_xy=xy[mtx.idx_b, :],
                   b_sizes=sizes[mtx.idx_b, :],
                   minimum_distance=minimum_distance)
    return mtx.fill(values=dist)


def distances(a_xy: NDArray, a_sizes: NDArray,
              b_xy: NDArray, b_sizes: NDArray) -> NDArray:
    """Shortest distance between rectangles

    Note: does not return negative distances for overlap. Overlap-> dist=0
    """

    assert a_xy.shape == a_sizes.shape and b_xy.shape == b_sizes.shape
    xy_dist = np.abs(a_xy - b_xy) - (a_sizes + b_sizes) / 2  # type: ignore
    xy_dist[np.where(xy_dist < 0)] = 0
    return np.hypot(xy_dist[:, 0], xy_dist[:, 0])


def distances_coordinate(rect_xy: NDArray, rect_sizes: NDArray,
                         coord_xy) -> NDArray:
    """distances between rectangles and coordinates

    Note: does not return negative distances for overlap. Overlap-> dist=0
    """
    return distances(a_xy=rect_xy, a_sizes=rect_sizes,
                     b_xy=coord_xy, b_sizes=np.array([0, 0]))


def distance_matrix(xy: NDArray, sizes: NDArray) -> NDArray:
    """Return matrix with distance between the rectangles"""
    assert xy.shape == sizes.shape

    mtx = CombinationMatrx(xy.shape[0])
    dist = distances(a_xy=xy[mtx.idx_a, :],
                     a_sizes=sizes[mtx.idx_a, :],
                     b_xy=xy[mtx.idx_b, :],
                     b_sizes=sizes[mtx.idx_b, :])
    return mtx.fill(values=dist)


def spatial_relation(a_xy: NDArray, a_sizes: NDArray,
                     b_xy: NDArray, b_sizes: Optional[NDArray]) -> NDArray:
    """spatial relation between two rectangles.
    returns angle and distance long the line of the spatial relation.
    """
    assert a_xy.shape == a_xy.shape
    rtn = cartesian2polar(b_xy - a_xy)   # type: ignore
    d_a = center_edge_distance(angles=rtn[:, 1], rect_sizes=a_sizes)
    if b_sizes is None:
        d_b = 0
    else:
        d_b = center_edge_distance(angles=rtn[:, 1], rect_sizes=b_sizes)
    rtn[:, 0] = rtn[:, 0] - d_a - d_b

    return rtn


def spatial_relation_coordinate(rect_xy: NDArray, rect_sizes: NDArray,
                                coord_xy: NDArray) -> NDArray:
    """spatial relation between rectangles and coordinate
    returns angle and distance long the line of the spatial relation.
    """
    return spatial_relation(a_xy=rect_xy, a_sizes=rect_sizes,
                            b_xy=coord_xy, b_sizes=None)


def spatial_relation_dot(rect_xy: NDArray, rect_sizes: NDArray,
                         dot_xy: NDArray, dot_diameter: NDArray) -> NDArray:
    """spatial relation between two rectangles.
    returns angle and distance long the line of the spatial relation.
    """
    assert (rect_xy.shape == rect_sizes.shape and
            dot_diameter.shape[0] == dot_xy.shape[0])
    # relation with dot center
    rtn = spatial_relation(a_xy=rect_xy, a_sizes=rect_sizes,
                           b_xy=dot_xy, b_sizes=None)
    rtn[:, 0] = rtn[:, 0] - dot_diameter  # type: ignore
    return rtn


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


def corner_distances(rect_xy: NDArray, rect_sizes: NDArray,
                     coord_xy: NDArray) -> NDArray:
    """ 2D array of euclidean distances (n, 4) of all corners of rectangles to
    coordinates"""

    assert rect_xy.shape == rect_sizes.shape == coord_xy.shape
    # set shapes for coord_xy (n, 2, 1)
    coord_xy = coord_xy.reshape(coord_xy.shape[0], 2, 1)
    xy_dist = corner_tensor(
        rect_xy=rect_xy, rect_sizes=rect_sizes) - coord_xy  # type: ignore
    return np.hypot(xy_dist[:, 0, :], xy_dist[:, 1, :])


# TODO CHECK this
def overlap_with_dots(rect_xy: NDArray, rect_sizes: NDArray,
                      dot_xy: NDArray, dot_diameter: NDArray,
                      minimum_distance: float = 0) -> Union[np.bool_, NDArray[np.bool_]]:
    """check if rectangles overlap with dots"""

    assert (rect_xy.shape == rect_sizes.shape == dot_xy.shape and
            dot_diameter.shape[0] == dot_xy.shape[0])

    dot_radii = dot_diameter.reshape((rect_xy.shape[0], 1)) / 2  # (n, 1)-array

    # overlap corner
    eucl_dist = corner_distances(rect_xy=rect_xy, rect_sizes=rect_sizes,
                                 coord_xy=dot_xy)
    if rect_xy.shape[0] > 1:
        overlap_corner = np.any(
            eucl_dist < dot_radii + minimum_distance, axis=1)
    else:
        overlap_corner = np.any(eucl_dist < dot_radii + minimum_distance)

    # overlap dot outer points
    # corners might not be in dot, but dot or intersects rect edge or is inside rect
    # check xy dist overlap
    dot_outer = np_dots.outer_points(dot_xy=dot_xy, dot_diameter=dot_diameter)
    rect_xy = rect_xy.reshape(rect_xy.shape[0], 2, 1)  # make 3D
    rect_sizes = rect_sizes.reshape(rect_sizes.shape[0], 2, 1)  # make 3D
    xy_diff = np.abs(dot_outer - rect_xy) - rect_sizes / 2  # type: ignore
    # check all true on xy of each outpoint, than  true  any of these comparison is true (for each point)
    comp = np.all(xy_diff < minimum_distance, axis=1)
    if comp.shape[0] > 1:
        overlap_outer = np.any(comp, axis=1)
    else:
        overlap_outer = np.any(comp)

    return overlap_corner | overlap_outer


def inside_rectangle(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates with rectangles are fully inside the rectangle/rectangles"""
    raise NotImplementedError()


def inside_dots(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates that rectangles are fully inside the dot/dots"""
    raise NotImplementedError()


def center_edge_distance(angles: NDArray, rect_sizes: NDArray) -> NDArray[np.floating]:
    """Distance between rectangle center and rectangle edge along the line
    in direction of `angle`.
    """

    if rect_sizes.ndim == 1:
        rect_sizes = np.ones((angles.shape[0], 1)) * rect_sizes

    l_inside = np.empty(len(angles))
    # find vertical relations
    v_rel = (np.pi-np.pi/4 >= abs(angles)) & (abs(angles) > np.pi/4)
    # vertical relation: in case line cut rectangle at the top or bottom corner
    i = np.flatnonzero(v_rel)
    l_inside[i] = rect_sizes[i, 1] / 2 * np.cos(np.pi/2 - angles[i])
    # horizontal relation: in case line cut rectangle at the left or right corner
    i = np.flatnonzero(~v_rel)
    l_inside[i] = rect_sizes[i, 0] / 2 * np.cos(angles[i])

    return l_inside


# GO ON HERE: TEST THESE FUNCTIONS


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

    assert(rects_a_xy.shape == rects_a_sizes.shape ==
           rects_b_xy.shape == rects_b_sizes.shape)

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

    assert(rects_xy.shape == rects_sizes.shape == dots_xy.shape and
           dots_xy.shape[0] == dots_diameter.shape[0])

    for arr in (rects_xy, rects_sizes, dots_xy):
        if arr.ndim > 1 and arr.shape[0] > 1:
            n_rows = arr.shape[0]
            break
    rects_xy = as_array2d_nrow(rects_xy, n_rows=n_rows)
    rects_sizes = as_array2d_nrow(rects_sizes, n_rows=n_rows)
    dots_xy = as_array2d_nrow(dots_xy, n_rows=n_rows)
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
