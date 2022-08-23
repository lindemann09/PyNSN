"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Union
import numpy as np
from numpy.typing import NDArray

from .np_tools import as_vector
from .geometry import corner_tensor, dots_outer_points, polar2cartesian
# all init-functions require 2D arrays


class CoordinateSpatRel(object):

    def __init__(self, a_xy: NDArray, b_xy: NDArray) -> None:
        if a_xy.ndim == 1:
            a_xy = a_xy.reshape((1, len(a_xy)))  # make 2d
        self._xy_diff = b_xy - a_xy  # type: ignore

    def spatial_relations(self):
        """distances and angles between coordinates"""
        # basically cartesian to polar
        return np.array([np.hypot(self._xy_diff[:, 0], self._xy_diff[:, 1]),
                        np.arctan2(self._xy_diff[:, 1], self._xy_diff[:, 0])]).T

    def distances(self) -> NDArray:
        """Euclidean distances between coordinates (a and b)

        Note: At least one xy parameter has to be a 2D array
        """
        return np.hypot(self._xy_diff[:, 0], self._xy_diff[:, 1])


class DotSpatRel(CoordinateSpatRel):

    def __init__(self, a_xy: NDArray, a_diameter: NDArray,
                 b_xy: NDArray, b_diameter: NDArray) -> None:

        if a_xy.ndim == 1:
            a_xy = a_xy.reshape((1, len(a_xy)))
        if b_xy.ndim == 1:
            b_xy = b_xy.reshape((1, len(b_xy)))

        assert a_xy.shape[0] == a_diameter.shape[0] and \
            b_xy.shape[0] == b_diameter.shape[0]

        super().__init__(a_xy=a_xy, b_xy=b_xy)
        self._radii_sum = (a_diameter + b_diameter) / 2

    def distances(self) -> NDArray:
        """Euclidean distances between dots (a and b)"""

        # distance between centers minus the radii
        return CoordinateSpatRel.distances(self) - self._radii_sum

    def spatial_relations(self) -> NDArray:
        """distances and angles between dots"""
        rtn = CoordinateSpatRel.spatial_relations(self)
        rtn[:, 0] = rtn[:, 0] - self._radii_sum
        return rtn

    def overlaps(self, minimum_distance: float = 0) -> Union[NDArray[np.bool_], np.bool_]:
        """True if dots overlap"""
        # columns both negative -> it overlaps
        comp = self.distances() < minimum_distance
        if comp.shape[0] > 1:
            return np.all(comp, axis=1)
        else:
            return np.all(comp)

    def required_displacements(self, minimum_distance: float = 0) -> NDArray:
        # FIXME not tested
        sr = self.spatial_relations()
        ok = sr[:, 0] >= minimum_distance
        sr[ok, :] = 0
        sr[~ok, 0] = minimum_distance - sr[~ok, 0]
        sr[~ok, 1] = sr[~ok, 1] + np.pi
        return sr


class DotCoordinateSpatRel(DotSpatRel):

    def __init__(self, dot_xy: NDArray, dot_diameter: NDArray,
                 coord_xy: NDArray) -> None:
        super().__init__(a_xy=dot_xy,
                         a_diameter=dot_diameter,
                         b_xy=coord_xy,
                         b_diameter=np.zeros(coord_xy.shape[0]))


class RectangleSpatRel(CoordinateSpatRel):

    def __init__(self, a_xy: NDArray, a_sizes: NDArray,
                 b_xy: NDArray, b_sizes: NDArray) -> None:
        assert a_sizes.shape == a_xy.shape and b_sizes.shape == b_xy.shape

        super().__init__(a_xy=a_xy, b_xy=b_xy)
        if a_sizes.ndim == 1:
            self.a_sizes = a_sizes.reshape((1, len(a_sizes)))  # make 2d
        else:
            self.a_sizes = a_sizes
        if b_sizes.ndim == 1:
            self.b_sizes = b_sizes.reshape((1, len(b_sizes)))  # make 2d
        else:
            self.b_sizes = b_sizes

        self._xy_diff_rect = np.abs(self._xy_diff) - \
            (self.a_sizes + self.b_sizes) / 2

    def overlaps(self, minimum_distance: float = 0) -> Union[NDArray[np.bool_], np.bool_]:
        """True if rectangles overlap"""
        # columns both negative -> it overlaps
        comp = self._xy_diff_rect < minimum_distance
        if comp.shape[0] > 1:
            return np.all(comp, axis=1)
        else:
            return np.all(comp)

    def distances(self) -> NDArray:
        """Shortest distance between rectangles

        Note: does not return negative distances for overlap. Overlap-> dist=0
        """

        xy_dist = np.copy(self._xy_diff_rect)
        xy_dist[np.where(xy_dist < 0)] = 0
        return np.hypot(xy_dist[:, 0], xy_dist[:, 1])

    def spatial_relations(self) -> NDArray:
        """spatial relation between two rectangles.
        returns distance and angle long the line of the spatial relation.
        """
        rtn = CoordinateSpatRel.spatial_relations(self)
        d_a = _center_edge_distance(angles=rtn[:, 1], rect_sizes=self.a_sizes)
        d_b = _center_edge_distance(angles=rtn[:, 1], rect_sizes=self.b_sizes)
        rtn[:, 0] = rtn[:, 0] - d_a - d_b

        return rtn

    def required_displacements(self, minimum_distance: float = 0) -> NDArray:
        """Minimum required displacement of rectangle B to have no overlap with
        rectangle A. Array with replacement vectors in polar coordinates.


        Returns:
            replacements as array of vectors in cartesian coordinates
            if replacements coordinate is np.nan no replacement required.

            Calculated movement direction (replacements[:, 0]) is always along the
            line between the two object center.
        """

        # spatial relations, but enlarge second rect by minimum distance
        sr = CoordinateSpatRel.spatial_relations(self)
        d_a = _center_edge_distance(angles=sr[:, 1],
                                    rect_sizes=self.a_sizes)
        d_b = _center_edge_distance(angles=sr[:, 1],
                                    rect_sizes=self.b_sizes + 2 * minimum_distance)
        sr[:, 0] = sr[:, 0] - d_a - d_b

        # set all nan and override overlapping relations
        rtn = np.full(sr.shape, np.nan)
        i = np.flatnonzero(sr[:, 0] < 0)
        sr[i, 1] = np.pi + sr[i, 1]  # opposite direction
        rtn[i, :] = polar2cartesian(sr[i, :])
        return rtn


class RectangleCoordinateSpatRel(RectangleSpatRel):

    def __init__(self, rect_xy: NDArray, rect_sizes: NDArray,
                 coord_xy: NDArray) -> None:
        super().__init__(a_xy=rect_xy,
                         a_sizes=rect_sizes,
                         b_xy=coord_xy,
                         b_sizes=np.zeros(coord_xy.shape))


class RectangleDotSpatRel(CoordinateSpatRel):

    def __init__(self, rect_xy: NDArray, rect_sizes: NDArray,
                 dot_xy: NDArray, dot_diameter: NDArray) -> None:

        if dot_xy.ndim == 1:
            dot_xy = dot_xy.reshape((1, len(dot_xy)))

        assert rect_xy.shape == rect_sizes.shape and \
            dot_xy.shape[0] == dot_diameter.shape[0]

        if rect_xy.ndim == 1:
            self.rect_xy = rect_xy.reshape((1, len(rect_xy)))  # make 2d
            self.rect_sizes = rect_sizes.reshape(
                (1, len(rect_sizes)))  # make 2d
        else:
            self.rect_xy = rect_xy
            self.rect_sizes = rect_sizes
        self.dot_xy = dot_xy

        super().__init__(a_xy=self.rect_xy, b_xy=self.dot_xy)
        # self.dot_diameter = dot_diameter
        self._dot_radii = as_vector(dot_diameter) / 2

    def spatial_relations(self) -> NDArray:
        """spatial relation between two rectangles.
        returns angle and distance long the line of the spatial relation.
        """

        # relation with dot center
        rtn = CoordinateSpatRel.spatial_relations(self)
        d_a = _center_edge_distance(angles=rtn[:, 1],
                                    rect_sizes=self.rect_sizes)
        rtn[:, 0] = rtn[:, 0] - d_a - self._dot_radii
        return rtn

    def distances(self) -> NDArray:
        return np.min(self._corner_distances(), axis=1)

    def _corner_distances(self) -> NDArray:
        # euclidean distance corner distances dot center shape=(n, 4)
        # set shapes for xy (n, 2, 1)
        dot_xy = self.dot_xy.reshape(self.dot_xy.shape[0], 2, 1)
        xy_dist = corner_tensor(rect_xy=self.rect_xy,
                                rect_sizes=self.rect_sizes) - dot_xy  # type: ignore
        return np.hypot(xy_dist[:, 0, :], xy_dist[:, 1, :])

    def overlaps(self, minimum_distance: float = 0) -> Union[np.bool_, NDArray[np.bool_]]:
        """check if rectangles overlap with dots"""

        dist_corner = self._corner_distances()
        # overlap corner
        if self.rect_xy.shape[0] > 1:
            # (n, 1)-array
            radii = self._dot_radii.reshape(self._dot_radii.shape[0], 1)
            overlap_corner = np.any(
                dist_corner < radii + minimum_distance, axis=1)
        else:
            overlap_corner = np.any(
                dist_corner < self._dot_radii + minimum_distance)

        # overlap dot outer points
        # corners might not be in dot, but dot or intersects rect edge or is inside rect
        # check xy dist overlap
        dot_outer = dots_outer_points(dot_xy=self.dot_xy,
                                      dot_radii=self._dot_radii)
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

    def required_displacements(self, minimum_distance: float = 0) -> NDArray:
        """Minimum required displacement of rectangle B to have no overlap with
        rectangle A. Array with replacement vectors in polar coordinates.

        Two objects do not overlap if replacements[:, 0] == 0.

        Returns:
            replacements as array of vectors in cartesian coordinates
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
        rtn = np.full(self._xy_diff.shape, np.nan)

        # only process overlapping items further
        min_spatial_rel = _min_spatial_relation(xy_diff=self._xy_diff[i_over, :],
                                                min_xy_dist=min_xy_dist[i_over, :])
        # subtract distance inside rectangles from euclidean distance
        l_inside_rect = _center_edge_distance(min_spatial_rel[:, 1],
                                              self.rect_sizes[i_over, :])
        min_spatial_rel[:, 0] = min_spatial_rel[:, 0] \
            - (l_inside_rect + self._dot_radii)
        # FIXME check if replacement is at all required here. spatial relation might be smaller than inside_rect + inside dot
        # insert calculated replacements into return array
        rtn[i_over, :] = polar2cartesian(min_spatial_rel)

        return rtn


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
    l_inside[i] = rect_sizes[i, 1] / (2 * np.cos(np.pi/2 - angles[i]))
    # horizontal relation: in case line cut rectangle at the left or right corner
    i = np.flatnonzero(~v_rel)
    l_inside[i] = rect_sizes[i, 0] / (2 * np.cos(angles[i]))
    return np.abs(l_inside)


def _min_spatial_relation(xy_diff: NDArray, min_xy_dist: NDArray) -> NDArray[np.floating]:
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
    x_smaller = np.abs(xy_diff[:, 0]) < np.abs(xy_diff[:, 1])
    xy_equal = xy_diff[:, 0] == xy_diff[:, 1]
    if np.sum(xy_equal) > 0:
        # randomly take some (via filtering) and add to x_smaller
        filter_vector = np.random.choice([False, True], size=len(xy_equal))
        x_smaller = x_smaller | (xy_equal & filter_vector)

    # calc required euclidean distances between object center
    min_eucl_dist = np.empty(xy_diff.shape[0])
    i = np.flatnonzero(x_smaller)  # x smaller -> set x to min_xy_dist
    min_eucl_dist[i] = np.hypot(min_xy_dist[i, 0], xy_diff[i, 1])
    i = np.flatnonzero(~x_smaller)  # x larger -> set y to min_xy_dist
    min_eucl_dist[i] = np.hypot(xy_diff[i, 0], min_xy_dist[i, 1])

    # return array of polar vectors
    # [[euclidean distances, directions between center]]
    directions = np.arctan2(xy_diff[:, 1], xy_diff[:, 0])
    return np.array([min_eucl_dist, directions]).T


def inside_rectangles(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates with dots are fully inside the rectangle/rectangles"""
    raise NotImplementedError()


def inside_dots(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates with dots are fully inside the dot/dots"""
    raise NotImplementedError()
