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


def line_point_othogonal(p1_line, p2_line, p3, outside_segment_nan=False):
    """project of P3 to line (p1, p2)
    Crossing point of the orthogonal to line (LP1, LP2) that goes to the
    point p3

    if outside_segment_nan, function checks if crosspoint is on the line segment
    (p1-p2) and if not values are set to nan.

    https://stackoverflow.com/questions/47177493/python-point-on-a-line-closest-to-third-point
    """
    p1 = np.atleast_2d(p1_line)
    dxy21 = np.atleast_2d(p2_line) - p1  # type: ignore
    dxy31 = np.atleast_2d(p3) - p1  # type: ignore
    a = np.sum(dxy21*dxy31, axis=1) / np.sum(dxy21**2, axis=1)
    cross_points = p1 + np.atleast_2d(a).T * dxy21
    if outside_segment_nan:
        # boolean array (n, 2=outside_x or y)
        outside = (np.maximum(p1_line, p2_line) < cross_points) |\
            (cross_points < np.minimum(p1_line, p2_line))
        outside = outside[:, 0] | outside[:, 1]  # vector
        cross_points[outside, :] = np.nan

    return cross_points


def corners(rect_xy: NDArray, rect_sizes: NDArray, lt_rb_only=False) -> NDArray:
    """tensor (n, 2, 4) with xy values of the four corners of the rectangles
    0=left-top, 1=right-top, 2=right-bottom, 3=left-bottom

    if lt_rb_only=True,
        return only  left-top and right-bottom point (n, 2_xy, 2=(lt, rb)
        """
    rect_xy = np.atleast_2d(rect_xy)
    rect_sizes2 = np.atleast_2d(rect_sizes) / 2
    right_top = rect_xy + rect_sizes2
    left_button = rect_xy - rect_sizes2

    if lt_rb_only:
        rtn = np.empty((rect_xy.shape[0], 2, 2))
        # left, top
        rtn[:, 0, 0] = left_button[:, 0]
        rtn[:, 1, 0] = right_top[:, 1]
        # right, bottom
        rtn[:, 0, 1] = right_top[:, 0]
        rtn[:, 1, 1] = left_button[:, 1]
    else:
        rtn = np.empty((rect_xy.shape[0], 2, 4))
        rtn[:, :, 1] = right_top  # right top
        rtn[:, :, 3] = left_button  # left bottom
        # left, top
        rtn[:, 0, 0] = left_button[:, 0]
        rtn[:, 1, 0] = right_top[:, 1]
        # right, bottom
        rtn[:, 0, 2] = right_top[:, 0]
        rtn[:, 1, 2] = left_button[:, 1]

    return rtn


def center_edge_distance(angles: NDArray, rect_sizes: NDArray) -> NDArray[np.floating]:
    """Distance between rectangle center and rectangle edge along the line
    in direction of `angle`.
    """
    l_inside = np.empty(len(angles))
    rect_sizes2 = rect_sizes / 2
    # find horizontal relations
    mod_angle = np.abs(angles) % np.pi  # scale 0 to PI
    i_h = (mod_angle > np.pi*.75) | (mod_angle < np.pi*.25)
    # horizontal relation: in case line cut rectangle at the left or right edge
    l_inside[i_h] = rect_sizes2[i_h, 0] / np.cos(mod_angle[i_h])
    # vertical relation: in case line cut rectangle at the top or bottom edge
    l_inside[~i_h] = rect_sizes2[~i_h, 1] / np.cos(np.pi/2 - mod_angle[~i_h])
    return np.abs(l_inside)


def distances_along_polar_radius(rho: NDArray, xy_distances: NDArray) -> NDArray:
    """Calculates the Euclidean distances along the polar radius (rho) that
    correspond to x and y distances along the cartesian x and y.

    rho: array of n angle
    xy_distance: array (n, 2) with x ([:,0]) and y ([:, 1]) distances

    Returns 2-D array with Euclidean distance
    """
    rtn = np.empty_like(xy_distances)
    # find the point on the line between center that correspond to the
    # target displacement distance at x or y axis
    # to get no overlap vertically:
    target_x_diff = np.abs(np.tan(np.pi/2 - rho)) * xy_distances[:, 0]
    rtn[:, 0] = np.hypot(target_x_diff, xy_distances[:, 1])
    # to get no overlap horizontal: target_y_diff
    target_y_diff = np.abs(np.tan(rho)) * xy_distances[:, 1]
    rtn[:, 1] = np.hypot(target_y_diff, xy_distances[:, 0])

    return rtn
