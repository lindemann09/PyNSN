"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Optional, Union
import numpy as np
from numpy.typing import NDArray

from .np_tools import CombinationMatrx, as_vector
from .np_coordinates import cartesian2polar
from . import np_dots

# all init-functions require 2D arrays


class RectangleRelations(object):

    def __init__(self, a_xy: NDArray, a_sizes: NDArray,
                 b_xy: NDArray, b_sizes: NDArray) -> None:
        assert a_xy.shape == a_xy.shape
        assert b_sizes.shape == b_xy.shape
        self._xy_diff = b_xy - a_xy  # type: ignore
        self.a_sizes = a_sizes
        self.b_sizes = b_sizes

    def overlap(self, minimum_distance: float = 0) -> Union[NDArray[np.bool_], np.bool_]:
        """True if rectangles overlap

        Note: At least one xy parameter has to be a 2D array
        """
        # columns both negative -> it overlaps
        xy_dist = np.abs(self._xy_diff) - (self.a_sizes + self.b_sizes) / 2
        comp = xy_dist < minimum_distance
        if comp.shape[0] > 1:
            return np.all(comp, axis=1)
        else:
            return np.all(comp)

    def distances(self) -> NDArray:
        """Shortest distance between rectangles

        Note: does not return negative distances for overlap. Overlap-> dist=0
        """

        xy_dist = np.abs(self._xy_diff) - (self.a_sizes + self.b_sizes) / 2
        xy_dist[np.where(xy_dist < 0)] = 0
        return np.hypot(xy_dist[:, 0], xy_dist[:, 0])

    def spatial_relation(self) -> NDArray:
        """spatial relation between two rectangles.
        returns angle and distance long the line of the spatial relation.
        """
        rtn = cartesian2polar(self._xy_diff)
        d_a = _center_edge_distance(angles=rtn[:, 1], rect_sizes=self.a_sizes)
        if isinstance(self.b_sizes, int):
            d_b = 0
        else:
            d_b = _center_edge_distance(
                angles=rtn[:, 1], rect_sizes=self.b_sizes)
        rtn[:, 0] = rtn[:, 0] - d_a - d_b

        return rtn

    def is_inside(self) -> NDArray:
        """bool array indicates with rectangles are fully inside the rectangle/rectangles"""
        raise NotImplementedError()

    def required_displacement(self, minimum_distance: float = 0) -> NDArray:
        """Minimum required displacement of rectangle B to have no overlap with
        rectangle A. Array with replacement vectors in polar coordinates.

        Two objects do not overlap if replacements[:, 0] == 0.

        Returns:
            replacements as array of vectors in polar coordinates
            if distance (replacements[:, 0]) is 0, no replacement required.

            Calculated movement direction (replacements[:, 0]) is always along the
            line between the two object center.
        """
        # TODO not tested

        #
        # find overlapping objects
        #
        min_xy_dist = ((self.a_sizes + self.b_sizes) / 2) + minimum_distance
        overlapping = np.all(np.abs(self._xy_diff) - min_xy_dist < 0, axis=1)
        i_over = np.flatnonzero(overlapping)

        #
        # calc dist and direction for overlapping objects
        #
        rtn = np.zeros(self._xy_diff.shape)

        # only process overlapping items further
        min_spatial_rel = _min_spatial_relation(xy_dist=self._xy_diff[i_over, :],
                                                min_xy_dist=min_xy_dist[i_over, :])
        # subtract distance inside rectangles from euclidean distance
        l_inside_a = _center_edge_distance(min_spatial_rel[:, 1],
                                           self.a_sizes[i_over, :])
        l_inside_b = _center_edge_distance(min_spatial_rel[:, 1],
                                           self.b_sizes[i_over, :])
        min_spatial_rel[:, 0] = min_spatial_rel[:, 0] - \
            (l_inside_a + l_inside_b)

        # insert calculated replacements into return array
        rtn[i_over, :] = min_spatial_rel

        return rtn


class RectangleCoordinateRelations(RectangleRelations):

    def __init__(self, rect_xy: NDArray, rect_sizes: NDArray,
                 coord_xy: NDArray) -> None:
        super().__init__(a_xy=rect_xy,
                         a_sizes=rect_sizes,
                         b_xy=coord_xy,
                         b_sizes=np.zeros(coord_xy.shape))

    def is_inside(self) -> NDArray:
        """bool array indicates with rectangles are fully inside the rectangle/rectangles"""
        raise RuntimeError("Rectangles can not be inside a points.")


class MatrixRectangleRelations(object):

    def __init__(self, xy: NDArray, sizes: NDArray):
        mtx = CombinationMatrx(xy.shape[0])
        self._rr = RectangleRelations(a_xy=xy[mtx.idx_a, :],
                                      a_sizes=sizes[mtx.idx_a, :],
                                      b_xy=xy[mtx.idx_b, :],
                                      b_sizes=sizes[mtx.idx_b, :])
        self._mtx = mtx

    def overlap(self, minimum_distance: float = 0) -> NDArray:
        """Matrix with overlaps (True/False) between the rectangles"""
        return self._mtx.fill(values=self._rr.overlap(minimum_distance))

    def distances(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self._mtx.fill(values=self._rr.distances())


class RectangleDotRelations(object):

    def __init__(self, rect_xy: NDArray, rect_sizes: NDArray,
                 dot_xy: NDArray, dot_diameter: NDArray) -> None:
        assert rect_xy.shape == rect_sizes.shape
        assert dot_xy.shape[0] == dot_diameter.shape[0]

        self.rect_xy = rect_xy
        self.dot_xy = dot_xy
        self._xy_diff = dot_xy - rect_xy   # type: ignore
        self.rect_sizes = rect_sizes
        # self.dot_diameter = dot_diameter
        self._dot_radii = as_vector(dot_diameter) / 2

    def spatial_relation_dot(self) -> NDArray:
        """spatial relation between two rectangles.
        returns angle and distance long the line of the spatial relation.
        """

        # relation with dot center
        rtn = cartesian2polar(self._xy_diff)   # type: ignore
        d_a = _center_edge_distance(
            angles=rtn[:, 1], rect_sizes=self.rect_sizes)
        rtn[:, 0] = rtn[:, 0] - d_a - self._dot_radii
        return rtn

    def overlap(self, minimum_distance: float = 0) -> Union[np.bool_, NDArray[np.bool_]]:
        """check if rectangles overlap with dots"""

        rc = RectangleCorners(rect_xy=self.rect_xy, rect_sizes=self.rect_sizes)
        # overlap corner
        eucl_dist = rc.distances(coord_xy=self.dot_xy)
        if self.rect_xy.shape[0] > 1:
            # (n, 1)-array
            radii = self._dot_radii.reshape(self._dot_radii.shape[0], 1)
            overlap_corner = np.any(
                eucl_dist < radii + minimum_distance, axis=1)
        else:
            overlap_corner = np.any(
                eucl_dist < self._dot_radii + minimum_distance)

        # overlap dot outer points
        # corners might not be in dot, but dot or intersects rect edge or is inside rect
        # check xy dist overlap
        dot_outer = np_dots.outer_points(
            dot_xy=self.dot_xy, dot_radii=self._dot_radii)
        rect_xy = self.rect_xy.reshape(self.rect_xy.shape[0], 2, 1)  # make 3D
        rect_sizes = self.rect_sizes.reshape(
            self.rect_sizes.shape[0], 2, 1)  # make 3D
        xy_diff = np.abs(dot_outer - rect_xy) - rect_sizes / 2  # type: ignore
        # check all true on xy of each outpoint, than  true  any of these comparison is true (for each point)
        comp = np.all(xy_diff < minimum_distance, axis=1)
        if comp.shape[0] > 1:
            overlap_outer = np.any(comp, axis=1)
        else:
            overlap_outer = np.any(comp)

        return overlap_corner | overlap_outer

    def is_inside(self) -> NDArray:
        """bool array indicates that rectangles are fully inside the dot/dots"""
        raise NotImplementedError()

    def required_displacement_with_dot(self, minimum_distance: float = 0) -> NDArray:
        """Minimum required displacement of rectangle B to have no overlap with
        rectangle A. Array with replacement vectors in polar coordinates.

        Two objects do not overlap if replacements[:, 0] == 0.

        Returns:
            replacements as array of vectors in polar coordinates
            if distance (replacements[:, 0]) is 0, no replacement required.

            Calculated movement direction (replacements[:, 0]) is always along the
            line between the two object center.
        """

        # FIXME not tested
        dots_rect_sizes = np.transpose(self._dot_radii *
                                       np.full(shape=(2, self._dot_radii.shape[0]),
                                               fill_value=2))

        #
        # find xy overlapping objects (simplified: treat dot as rect)
        #
        min_xy_dist = ((self.rect_sizes + dots_rect_sizes) /
                       2) + minimum_distance
        overlapping = np.all(np.abs(self._xy_diff) - min_xy_dist < 0, axis=1)
        i_over = np.flatnonzero(overlapping)

        #
        # calc dist and direction for overlapping objects
        #
        rtn = np.zeros(self._xy_diff.shape)

        # only process overlapping items further
        min_spatial_rel = _min_spatial_relation(xy_dist=self._xy_diff[i_over, :],
                                                min_xy_dist=min_xy_dist[i_over, :])
        # subtract distance inside rectangles from euclidean distance
        l_inside_rect = _center_edge_distance(min_spatial_rel[:, 1],
                                              self.rect_sizes[i_over, :])
        min_spatial_rel[:, 0] = min_spatial_rel[:, 0] \
            - (l_inside_rect + self._dot_radii)
        # FIXME check if replacement is at all required here. spatial relation might be smaller than inside_rect + inside dot
        # insert calculated replacements into return array
        rtn[i_over, :] = min_spatial_rel

        return rtn


class RectangleCorners(object):

    def __init__(self, rect_xy: NDArray, rect_sizes: NDArray) -> None:
        """Corners is a tensor (n, 2, 4) with xy values of the four corners"""
        corner = np.empty((rect_xy.shape[0], 2, 4))

        rect_sizes2 = rect_sizes / 2
        corner[:, :, 1] = rect_xy + rect_sizes2  # right top
        corner[:, :, 3] = rect_xy - rect_sizes2  # left bottom
        # left, top
        corner[:, 0, 0] = corner[:, 0, 3]
        corner[:, 1, 0] = corner[:, 1, 1]
        # right, bottom
        corner[:, 0, 2] = corner[:, 0, 1]
        corner[:, 1, 2] = corner[:, 1, 3]

        self.corner = corner

    def distances(self, coord_xy: NDArray) -> NDArray:
        """ 2D array of euclidean distances (n, 4) of all corners of the
        rectangles to coordinates"""

        # set shapes for coord_xy (n, 2, 1)
        coord_xy = coord_xy.reshape(coord_xy.shape[0], 2, 1)
        xy_dist = self.corner - coord_xy  # type: ignore
        return np.hypot(xy_dist[:, 0, :], xy_dist[:, 1, :])


def _center_edge_distance(angles: NDArray, rect_sizes: NDArray) -> NDArray[np.floating]:
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
