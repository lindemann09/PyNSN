"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import NDArray, ArrayLike

from . import rng

# all functions are 2D arrays (at least) as fist arguments

NORTH = np.pi/2
SOUTH = np.pi/-2
NORTH_WEST = np.pi*.75
NORTH_EAST = np.pi*.25


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


def corners(rect_xy: NDArray, rect_sizes_div2: NDArray, lt_rb_only=False) -> NDArray:
    """tensor (n, 2, 4) with xy values of the four corners of the rectangles
    0=left-top, 1=right-top, 2=right-bottom, 3=left-bottom

    if lt_rb_only=True,
        return only  left-top and right-bottom point (n, 2_xy, 2=(lt, rb)
        """
    rect_xy = np.atleast_2d(rect_xy)
    rect_sizes_div2 = np.atleast_2d(rect_sizes_div2)
    right_top = rect_xy + rect_sizes_div2
    left_button = rect_xy - rect_sizes_div2  # type: ignore

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


def center_edge_distance(angles: NDArray, rect_sizes_div2: NDArray) -> NDArray[np.floating]:
    """Distance between rectangle center and rectangle edge along the line
    in direction of `angle`.
    """
    l_inside = np.empty(len(angles))
    # find horizontal relations
    mod_angle = np.abs(angles) % np.pi  # scale 0 to PI
    i_h = (mod_angle > NORTH_WEST) | (mod_angle < NORTH_EAST)
    # horizontal relation: in case line cut rectangle at the left or right edge
    l_inside[i_h] = rect_sizes_div2[i_h, 0] / np.cos(mod_angle[i_h])
    # vertical relation: in case line cut rectangle at the top or bottom edge
    l_inside[~i_h] = rect_sizes_div2[~i_h, 1] / np.cos(NORTH - mod_angle[~i_h])
    return np.abs(l_inside)


def chord_length(r: NDArray[np.floating], d: NDArray[np.floating])-> NDArray[np.floating]:
    """calculate the length of chord in a circle
    r: radius
    d: distance to center"""
    return 2.0 * np.sqrt(r**2 - d**2)

def chord_distance_to_center(r: NDArray, chord_len: NDArray)-> NDArray[np.floating]:
    """calculate distance to center of chord with a particular length
    r: radius
    d: distance to center"""
    return np.sqrt(r**2 - 0.25 * chord_len**2)