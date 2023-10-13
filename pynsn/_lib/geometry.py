"""geometry function"""

__author__ = "Oliver Lindemann <lindemann@cognitive-psychology.eu>"

from typing import List, Tuple, Union
import math

Coord2D = Tuple[float, float]
Coord2DLike = Union[Coord2D, List[float]]


def polar2cartesian(polar: Coord2DLike) -> Coord2D:
    """cartesian coordinates (x,y)
    polar has to be tuple representing polar coordinates (radius, angle)"""
    return polar[0] * math.cos(polar[1]), polar[0] * math.sin(polar[1])


def cartesian2polar(xy: Coord2DLike) -> Coord2D:
    """polar coordinates (radius, angle)
    xy has to be a tuple of two numerals
    """
    return math.hypot(xy[0], xy[1]), math.atan2(xy[1], xy[0])


def cartesian2image_coordinates(xy: Coord2DLike, image_size: Coord2DLike) -> Coord2D:
    """convert 2D cartesian to image coordinates with (0,0) at top left and
    reversed y axis
    """
    return xy[0] + image_size[0] / 2, xy[1] * -1 + image_size[1] / 2


# def jitter_identical_coordinates(xy: NDArray, jitter_size: float = 0.1) -> NDArray:
#     """jitters points with identical position"""

#     for idx, ref_object in enumerate(xy):
#         # find identical positions
#         identical = np.flatnonzero(np.all(np.equal(xy, ref_object), axis=1))
#         if len(identical) > 1:
#             for x in identical:  # jitter all identical positions
#                 if x != idx:
#                     tmp = np.atleast_2d((jitter_size,
#                                         rng.generator.random() * 2 * np.pi))
#                     xy[x, :] = xy[x, :] - polar2cartesian(tmp)[0]

#     return xy


# def corners(rect_xy: NDArray, rect_sizes_div2: NDArray, lt_rb_only=False) -> NDArray:
#     """tensor (n, 2, 4) with xy values of the four corners of the rectangles
#     0=left-top, 1=right-top, 2=right-bottom, 3=left-bottom

#     if lt_rb_only=True,
#         return only  left-top and right-bottom point (n, 2_xy, 2=[lt, rb])
#         """
#     rect_xy = np.atleast_2d(rect_xy)
#     rect_sizes_div2 = np.atleast_2d(rect_sizes_div2)
#     right_top = rect_xy + rect_sizes_div2
#     left_button = rect_xy - rect_sizes_div2  # type: ignore

#     if lt_rb_only:
#         rtn = np.empty((rect_xy.shape[0], 2, 2))
#         # left, top
#         rtn[:, 0, 0] = left_button[:, 0]
#         rtn[:, 1, 0] = right_top[:, 1]
#         # right, bottom
#         rtn[:, 0, 1] = right_top[:, 0]
#         rtn[:, 1, 1] = left_button[:, 1]
#     else:
#         rtn = np.empty((rect_xy.shape[0], 2, 4))
#         rtn[:, :, 1] = right_top  # right top
#         rtn[:, :, 3] = left_button  # left bottom
#         # left, top
#         rtn[:, 0, 0] = left_button[:, 0]
#         rtn[:, 1, 0] = right_top[:, 1]
#         # right, bottom
#         rtn[:, 0, 2] = right_top[:, 0]
#         rtn[:, 1, 2] = left_button[:, 1]

#     return rtn


# def center_edge_distance(angles: NDArray, rect_sizes_div2: NDArray) -> NDArray[np.floating]:
#     """Distance between rectangle center and rectangle edge along the line
#     in direction of `angle`.
#     """
#     l_inside = np.empty(len(angles))
#     # find horizontal relations
#     mod_angle = np.abs(angles) % np.pi  # scale 0 to PI
#     i_h = (mod_angle > NORTH_WEST) | (mod_angle < NORTH_EAST)
#     # horizontal relation: in case line cut rectangle at the left or right edge
#     l_inside[i_h] = rect_sizes_div2[i_h, 0] / np.cos(mod_angle[i_h])
#     # vertical relation: in case line cut rectangle at the top or bottom edge
#     l_inside[~i_h] = rect_sizes_div2[~i_h, 1] / np.cos(NORTH - mod_angle[~i_h])
#     return np.abs(l_inside)


# def chord_length(r: NDArray[np.floating], d: NDArray[np.floating]) -> NDArray[np.floating]:
#     """calculate the length of chord in a circle
#     r: radius
#     d: distance to center"""
#     return 2.0 * np.sqrt(r**2 - d**2)


# def chord_distance_to_center(r: NDArray, chord_len: NDArray) -> NDArray[np.floating]:
#     """calculate distance to center of chord with a particular length
#     r: radius
#     d: distance to center"""
#     return np.sqrt(r**2 - 0.25 * chord_len**2)


# def line_circle_intersection(line_points: NDArray[np.floating],
#                              line_directions: NDArray[np.floating],
#                              circle_center: NDArray[np.floating],
#                              circle_radii: NDArray,
#                              both_intersections: bool = False) -> NDArray[np.floating]:
#     """
#     Intersections of lines with the circles in one or both directions of the
#     line.

#     Parameters
#     ----------
#     line_points : 2d-array
#     line_rho : 1d-array of directions
#     circle_center : 2d-array
#     circle_radii: radius of circle
#     both_intersections: boolean (default=False)
#         if True, it returns both intersection points, otherwise only the
#         intersection in line direction

#     Returns
#     -------
#     intersecting points: NDArray
#         if both_intersections == True:
#             NDArray (n, 4)
#             The two points [x_1,y_1, x_2, y_2] per row

#         if both_intersections == False:
#             NDArray (n, 2)
#             only intersection points in the direction of rho

#     References
#     ----------
#     http://mathworld.wolfram.com/Circle-LineIntersection.html

#     https://github.com/ajhynes7/scikit-spatial/blob/master/src/skspatial/objects/circle.py


#     """

#     circle_center = np.atleast_2d(circle_center)
#     points_a = np.atleast_2d(line_points) - circle_center
#     points_b = points_a + polar2cartesian(
#         np.column_stack((np.ones_like(line_directions), line_directions)))

#     d = points_b - points_a

#     # Pre-compute variables common to x and y equations.
#     d_r_squared = d[:, 0]**2 + d[:, 1]**2
#     determinant = points_a[:, 0] * points_b[:, 1] - \
#         points_a[:, 1] * points_b[:, 0]
#     discriminant = np.squeeze(circle_radii)**2 * d_r_squared - determinant**2
#     discriminant[discriminant < 0] = np.nan

#     root = np.sqrt(discriminant)
#     x_sign = np.where(d[:, 1] < 0, -1, 1)  # sign with zero -> 1

#     if both_intersections:
#         points_ab = np.column_stack((
#             determinant * d[:, 1] - x_sign * d[:, 0] * root,  # x1
#             -determinant * d[:, 0] - np.abs(d[:, 1]) * root,  # y1
#             determinant * d[:, 1] + x_sign * d[:, 0] * root,  # x2
#             -determinant * d[:, 0] + np.abs(d[:, 1]) * root  # y2
#         ))
#         circle_center = np.column_stack((circle_center, circle_center))
#     else:
#         sign = np.where(np.sin(line_directions) < 0, -1, 1)  # which point
#         points_ab = np.column_stack((
#             determinant * d[:, 1] + sign * x_sign * d[:, 0] * root,
#             -determinant * d[:, 0] + sign * np.abs(d[:, 1]) * root
#         ))

#     return (points_ab / np.atleast_2d(d_r_squared).T) + circle_center
