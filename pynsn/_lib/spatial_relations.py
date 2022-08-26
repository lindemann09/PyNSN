"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Union
import numpy as np
from numpy.typing import NDArray

from .np_tools import as_vector, index_minmum_rowwise
from .geometry import polar2cartesian
from .._lib import geometry
# all init-functions require 2D arrays


class CoordinateSpatRel(object):

    def __init__(self, a_xy: NDArray, b_xy: NDArray) -> None:
        self._xy_diff = np.atleast_2d(b_xy) \
            - np.atleast_2d(a_xy)  # type: ignore

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

        assert a_xy.shape == b_xy.shape and \
            a_xy.shape[0] == a_diameter.shape[0] == b_diameter.shape[0]

        super().__init__(a_xy=a_xy, b_xy=b_xy)
        self._radii_sum = (as_vector(a_diameter) + as_vector(b_diameter)) / 2

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
        assert a_sizes.shape == a_xy.shape == b_sizes.shape == b_xy.shape

        super().__init__(a_xy=a_xy, b_xy=b_xy)
        self.a_sizes = np.atleast_2d(a_sizes)
        self.b_sizes = np.atleast_2d(b_sizes)
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
        self.rect_xy = np.atleast_2d(rect_xy)
        self.rect_sizes = np.atleast_2d(rect_sizes)
        self.dot_xy = np.atleast_2d(dot_xy)

        self._dot_radii = as_vector(dot_diameter) / 2
        self.__corners = None
        self.__corner_rel = None
        self.__edge_rel = None

        assert self.rect_xy.shape == self.rect_sizes.shape == self.dot_xy.shape and \
            self.dot_xy.shape[0] == self._dot_radii.shape[0]

        super().__init__(a_xy=self.rect_xy, b_xy=self.dot_xy)

    @property
    def _corners(self) -> NDArray:
        """tensor (n, 2, 4) with xy values of the four corners of the rectangles
        0=left-top, 1=right-top, 2=right-bottom, 3=left-bottom
        """
        if self.__corners is None:
            self.__corners = np.empty((self.rect_xy.shape[0], 2, 4))
            rect_sizes2 = self.rect_sizes / 2
            self.__corners[:, :, 1] = self.rect_xy + rect_sizes2  # right top
            self.__corners[:, :, 3] = self.rect_xy - rect_sizes2  # left bottom
            # left, top
            self.__corners[:, 0, 0] = self.__corners[:, 0, 3]
            self.__corners[:, 1, 0] = self.__corners[:, 1, 1]
            # right, bottom
            self.__corners[:, 0, 2] = self.__corners[:, 0, 1]
            self.__corners[:, 1, 2] = self.__corners[:, 1, 3]

        return self.__corners

    @property
    def _corner_relations(self) -> NDArray:
        """Euclidean distance and angle of nearest corner and dot centers
        shape=(n, 2)
        """

        if self.__corner_rel is None:
            xy_dist = self._corners - \
                np.atleast_3d(self.dot_xy)  # type: ignore
            distances = np.hypot(xy_dist[:, 0, :], xy_dist[:, 1, :])
            # find nearest corner per rect dist.shape=(n, 4)
            # Problem: one rect could have multiple corners with min distance
            #   --> choose randomly just one near corner of each rect
            idx_nc = index_minmum_rowwise(distances)

            # spatial relations and nearest corners
            # return array contains only one nearest corner per rect (n, 2=parameter)
            distances = distances[idx_nc[:, 0], idx_nc[:, 1]] - self._dot_radii
            directions = np.arctan2(
                xy_dist[idx_nc[:, 0], 1, idx_nc[:, 1]],
                xy_dist[idx_nc[:, 0], 0, idx_nc[:, 1]])

            self.__corner_rel = np.array([distances, directions]).T

        return self.__corner_rel

    @property
    def _edge_relations(self) -> NDArray:
        """spatial relations between nearest edge and dot center"""

        if self.__edge_rel is None:
            edge_points = np.array([(0, 1),  # north
                                    (1, 2),  # east
                                    (3, 2),  # south
                                    (0, 3)  # west
                                    ])
            # find nearest edge cross points to dot centers
            # (ecp -> project dot center to rect edge)
            # ecp_xy_diff (n, 2=xy, 4=edges)
            ecp_xy_diff = np.empty_like(self._corners)
            for i in range(4):
                # diff(edge, dot center): project of dot center to edge
                ecp_xy_diff[:, :, i] = geometry.line_point_othogonal(
                    p1_line=self._corners[:, :, edge_points[i, 0]],
                    p2_line=self._corners[:, :, edge_points[i, 1]],
                    p3=self.dot_xy) - self.dot_xy

            # Euclidian distance of each edge with dot center (n, 4=edges)
            distances = np.hypot(ecp_xy_diff[:, 0, :],
                                 ecp_xy_diff[:, 1, :])
            # find nearest edge
            idx_ne = index_minmum_rowwise(distances)
            # return array contains nearest edge per rect (n, 2=parameter)
            distances = distances[idx_ne[:, 0], idx_ne[:, 1]] - self._dot_radii
            directions = np.arctan2(
                ecp_xy_diff[idx_ne[:, 0], 1, idx_ne[:, 1]],
                ecp_xy_diff[idx_ne[:, 0], 0, idx_ne[:, 1]])

            self.__edge_rel = np.array([distances, directions]).T

        return self.__edge_rel

    def spatial_relations(self) -> NDArray:
        """spatial relation between rectangles and dots.
        returns angle and distance long the line of the spatial relation.
        """
        # FIXME not tested.
        # find closer corner idx
        idx_cc = self._corner_relations[:, 0] < self._edge_relations[:, 0]
        # spat rel shape =(n, 2=parameter)
        spat_rel = np.empty((self.rect_xy.shape[0], 2))
        spat_rel[idx_cc, :] = self._corner_relations[idx_cc, :]
        spat_rel[~idx_cc, :] = self._edge_relations[~idx_cc, :]
        return spat_rel

    def distances(self) -> NDArray:
        # FIXME not crrect does not take intp account edbge relation
        return np.min(self._corner_relations[:, :, 0], axis=1) \
            - self._dot_radii

    def is_corner_inside(self, minimum_distance: float = 0) -> Union[np.bool_, NDArray[np.bool_]]:
        """check if rectangles corner overlap with dots"""
        # overlap corner
        if self.rect_xy.shape[0] > 1:
            # (n, 1)-array
            return np.any(
                self._corner_relations[:, :, 0] < np.atleast_2d(self._dot_radii).T + minimum_distance, axis=1)
        else:
            return np.any(
                self._corner_relations[:, :, 0] < self._dot_radii + minimum_distance)

    def is_edge_intersecting(self, minimum_distance: float = 0) -> Union[np.bool_, NDArray[np.bool_]]:
        """check if rectangles edge intersect with dots"""
        # --> overlap dot outer points
        # corners might not be in dot, but dot or intersects rect edge or is inside rect
        # check xy dist overlap
        # FIXME USED self._edge_realtions!
        xy_diff = self._edge_relations(xy_distance_only=True)
        # check all true on xy of each outpoint, than  true  any of these comparison is true (for each point)
        comp = np.all(xy_diff < minimum_distance, axis=1)
        if comp.shape[0] > 1:
            return np.any(comp, axis=1)
        else:
            return np.any(comp)

    def overlaps(self, minimum_distance: float = 0) -> Union[np.bool_, NDArray[np.bool_]]:
        """check if rectangles overlap with dots"""

        return self.is_corner_inside() | self.is_edge_intersecting()

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
        # incorrect
        sr = self.spatial_relations()
        sr[:, 0] = sr[:, 0] - minimum_distance

        # set all nan and override overlapping relations
        rtn = np.full(sr.shape, np.nan)
        i = np.flatnonzero(sr[:, 0] < 0)
        sr[i, 1] = np.pi + sr[i, 1]  # opposite direction
        rtn[i, :] = polar2cartesian(sr[i, :])
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


# def _min_spatial_relation(xy_diff: NDArray, min_xy_dist: NDArray) -> NDArray[np.floating]:
#     """Calculate required minimum Euclidean distance and direction
#     (polar coordinates) from xy distances and min_minimum required xy distances

#     dist(x)**2 + dist(y)**2 = euclidean_dist(xy)**2
#     set e.g. x to minimum dist:
#       sqrt(min_dist(x)**2 + dist(y)**2) = min_euclidean_dist(xy)

#     Set the smaller xy distance to minimum dist and calc required movement
#     long the line between object centers (Euclidean distances minus length of
#     line inside objects).
#     """

#     # find x smaller then y distance
#     # if equals, randomly choose some and treat them like x_smaller
#     x_smaller = np.abs(xy_diff[:, 0]) < np.abs(xy_diff[:, 1])
#     xy_equal = xy_diff[:, 0] == xy_diff[:, 1]
#     if np.sum(xy_equal) > 0:
#         # randomly take some (via filtering) and add to x_smaller
#         filter_vector = np.random.choice([False, True], size=len(xy_equal))
#         x_smaller = x_smaller | (xy_equal & filter_vector)

#     # calc required euclidean distances between object center
#     min_eucl_dist = np.empty(xy_diff.shape[0])
#     i = np.flatnonzero(x_smaller)  # x smaller -> set x to min_xy_dist
#     min_eucl_dist[i] = np.hypot(min_xy_dist[i, 0], xy_diff[i, 1])
#     i = np.flatnonzero(~x_smaller)  # x larger -> set y to min_xy_dist
#     min_eucl_dist[i] = np.hypot(xy_diff[i, 0], min_xy_dist[i, 1])

#     # return array of polar vectors
#     # [[euclidean distances, directions between center]]
#     directions = np.arctan2(xy_diff[:, 1], xy_diff[:, 0])
#     return np.array([min_eucl_dist, directions]).T


def inside_rectangles(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates with dots are fully inside the rectangle/rectangles"""
    raise NotImplementedError()


def inside_dots(xy: NDArray, sizes: NDArray) -> NDArray:
    """bool array indicates with dots are fully inside the dot/dots"""
    raise NotImplementedError()


# def dots_cardinal_points(dot_xy: NDArray, dot_radii: NDArray) -> NDArray:
#     """returns a tensor (n, 2, 4) with four outer points (left, top, right, bottom)
#     on the cardinal axes"""
#     assert dot_xy.shape[0] == dot_radii.shape[0]

#     # make tensor (n, 2, 4) with dot_xy at each [:, :, x]
#     rtn = dot_xy.reshape((dot_xy.shape[0], dot_xy.shape[1], 1)) * \
#         np.ones((dot_xy.shape[0], 2, 4))
#     rtn[:, 0, 0] = rtn[:, 0, 0] - dot_radii  # left
#     rtn[:, 1, 1] = rtn[:, 1, 1] + dot_radii  # top
#     rtn[:, 0, 2] = rtn[:, 0, 2] + dot_radii  # right
#     rtn[:, 1, 3] = rtn[:, 1, 3] - dot_radii  # bottom
#     return rtn
