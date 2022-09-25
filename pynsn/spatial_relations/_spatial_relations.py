"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod

import warnings
from typing import Union
import numpy as np
from numpy.typing import NDArray

from .._lib.np_tools import find_value_rowwise
from .._lib.geometry import polar2cartesian
from .._lib import geometry

from .._object_arrays.dot_array import BaseDotArray
from .._object_arrays.rectangle_array import BaseRectangleArray
# all init-functions require 2D arrays


class ABCSpatialRelations(metaclass=ABCMeta):

    @abstractmethod
    def overlaps(self):
        """True is objects overlap"""

# FIXME TODO
#    @abstractmethod
#    def is_inside(self):
#        """True is objects are inside"""

    @abstractmethod
    def spatial_relations(self):
        """Distances and angles between objects"""

    @abstractmethod
    def required_displacements(self):
        """Distances and angles between objects"""


class CoordinateCoordinate(object):

    def __init__(self,
                 a_xy: NDArray,
                 b_xy: NDArray) -> None:
        self._xy_diff = np.atleast_2d(b_xy) \
            - np.atleast_2d(a_xy)  # type: ignore
        self._spatial_relations = None

    def spatial_relations(self):
        """distances and angles between coordinates"""
        # basically cartesian to polar
        if self._spatial_relations is None:
            self._spatial_relations = np.array([
                np.hypot(self._xy_diff[:, 0], self._xy_diff[:, 1]),
                np.arctan2(self._xy_diff[:, 1], self._xy_diff[:, 0])]).T
        return self._spatial_relations


class DotDot(CoordinateCoordinate, ABCSpatialRelations):

    def __init__(self,
                 dots_a: BaseDotArray,
                 dots_b: BaseDotArray) -> None:
        super().__init__(a_xy=dots_a.xy, b_xy=dots_b.xy)
        self._radii_sum = (dots_a.diameter + dots_b.diameter) / 2
        assert self._xy_diff.shape[0] == self._radii_sum.shape[0]

    def distances(self) -> NDArray:
        """Euclidean distances between dots (a and b)"""

        # distance between centers minus the radii
        return self.spatial_relations()[:, 0]

    def spatial_relations(self) -> NDArray:
        """distances and angles between dots (cached return data)"""
        if self._spatial_relations is None:
            rel = CoordinateCoordinate.spatial_relations(self)
            rel[:, 0] = rel[:, 0] - self._radii_sum
            self._spatial_relations = rel
        return self._spatial_relations

    def overlaps(self, minimum_distance: float = 0) -> Union[NDArray[np.bool_], np.bool_]:
        """True if dots overlap"""
        # columns both negative -> it overlaps
        comp = self.spatial_relations()[:, 0] < minimum_distance
        if comp.shape[0] > 1:
            return np.all(comp, axis=1)
        else:
            return np.all(comp)

    def required_displacements(self, minimum_distance: float = 0) -> NDArray:
        """Minimum required displacement of rectangle B to have no overlap with
        rectangle A.

        Two objects do not overlap if replacements[:, 0] == 0.

        Returns:
            replacements as array of vectors in cartesian coordinates
            if distance (replacements[:, 0]) is 0, no replacement required.

            Calculated movement direction (replacements[:, 0]) is always along the
            line between the two object center.
        """
        # FIXME not tested
        sr = self.spatial_relations()
        sr[:, 0] = sr[:, 0] - minimum_distance

        i = np.flatnonzero(sr[:, 0] < 0)  # overlaps
        sr[i, 1] = np.pi + sr[i, 1]  # opposite direction
        # set all nan and override overlapping relations
        rtn = np.full(sr.shape, np.nan)
        rtn[i, :] = polar2cartesian(sr[i, :])
        return rtn


class DotCoordinate(DotDot, ABCSpatialRelations):

    def __init__(self,
                 dots: BaseDotArray,
                 coord_xy: NDArray) -> None:
        dots_b = BaseDotArray(xy=coord_xy,
                              diameter=np.zeros(coord_xy.shape[0]))
        super().__init__(dots_a=dots, dots_b=dots_b)


class RectangleRectangle(CoordinateCoordinate, ABCSpatialRelations):

    def __init__(self,
                 rectangles_a: BaseRectangleArray,
                 rectangles_b: BaseRectangleArray) -> None:

        super().__init__(a_xy=rectangles_a.xy, b_xy=rectangles_b.xy)
        self.a_sizes = rectangles_a.sizes
        if np.all(rectangles_b.sizes == 0):  # Coordinate
            self.b_sizes = None
            self._xy_diff_rect = np.abs(self._xy_diff) - self.a_sizes / 2
        else:
            self.b_sizes = rectangles_b.sizes
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

    def xy_distances(self) -> NDArray:
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
        if self._spatial_relations is None:
            rel = CoordinateCoordinate.spatial_relations(self)
            d_a = _center_edge_distance(angles=rel[:, 1],
                                        rect_sizes=self.a_sizes)
            if self.b_sizes is None:
                d_b = 0
            else:
                d_b = _center_edge_distance(angles=rel[:, 1],
                                            rect_sizes=self.b_sizes)
            rel[:, 0] = rel[:, 0] - d_a - d_b
            self._spatial_relations = rel
        return self._spatial_relations

    def required_displacements(self, minimum_distance: float = 0) -> NDArray:
        """Minimum required displacement of rectangle B to have no overlap with
        rectangle A.

        Returns:
            replacements as array of vectors in cartesian coordinates
            if replacements coordinate is np.nan no replacement required.

            Calculated movement direction (replacements[:, 0]) is always along the
            line between the two object center.
        """

        # spatial relations, but enlarge first rect by minimum distance
        sr = CoordinateCoordinate.spatial_relations(self)
        d_a = _center_edge_distance(angles=sr[:, 1],
                                    rect_sizes=self.a_sizes + 2 * minimum_distance)
        if self.b_sizes is None:
            d_b = 0
        else:
            d_b = _center_edge_distance(angles=sr[:, 1],
                                        rect_sizes=self.b_sizes)
        sr[:, 0] = sr[:, 0] - d_a - d_b

        # set all nan and override overlapping relations
        rtn = np.full(sr.shape, np.nan)
        i = np.flatnonzero(sr[:, 0] < 0)
        sr[i, 1] = np.pi + sr[i, 1]  # opposite direction
        rtn[i, :] = polar2cartesian(sr[i, :])
        return rtn


class RectangleCoordinate(RectangleRectangle, ABCSpatialRelations):

    def __init__(self,
                 rectangles: BaseRectangleArray,
                 coord_xy: NDArray) -> None:
        rects_b = BaseRectangleArray(xy=coord_xy,
                                     sizes=np.zeros(coord_xy.shape))
        super().__init__(rectangles_a=rectangles, rectangles_b=rects_b)


class RectangleDot(CoordinateCoordinate, ABCSpatialRelations):

    def __init__(self,
                 rectangles: BaseRectangleArray,
                 dots: BaseDotArray) -> None:

        self.rect_xy = rectangles.xy
        self.rect_sizes = rectangles.sizes
        self.dot_xy = dots.xy
        self.dot_radii = dots.diameter / 2

        n_rects = len(self.rect_xy)
        n_dots = len(self.dot_xy)
        if n_rects > 1 and n_dots == 1:
            self.dot_xy = self.dot_xy * np.ones((n_rects, 1))
            self.dot_radii = self.dot_radii * np.ones(n_rects)
        elif n_rects == 1 and n_dots > 1:
            ones = np.ones((n_dots, 1))
            self.rect_xy = self.rect_xy * ones
            self.rect_sizes = self.rect_sizes * ones

        assert self.rect_xy.shape == self.dot_xy.shape

        self.__corners = None
        self._spatial_relations = None

        super().__init__(a_xy=self.rect_xy, b_xy=self.dot_xy)

    @property
    def _corners(self) -> NDArray:
        """tensor (n, 2, 4) with xy values of the four corners of the rectangles
        0=left-top, 1=right-top, 2=right-bottom, 3=left-bottom
        """
        if self.__corners is None:
            self.__corners = geometry.corners(rect_xy=self.rect_xy,
                                              rect_sizes=self.rect_sizes,
                                              lt_rb_only=False)
        return self.__corners

    def corner_relations(self, nearest_corners: bool = True) -> NDArray:
        """Euclidean distance and angle of nearest or farthest corners and dot centers
        shape=(n, 2)
        """

        xy_dist = self._corners - np.atleast_3d(self.dot_xy)  # type: ignore

        distances = np.hypot(xy_dist[:, 0, :], xy_dist[:, 1, :])
        # find nearest corner per rect dist.shape=(n, 4)
        # Problem: one rect could have multiple corners with min distance
        #   --> choose randomly just one near corner of each rect
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', 'All-NaN slice encountered')
            if nearest_corners:
                values = np.nanmin(distances, axis=1)
            else:
                values = np.nanmax(distances, axis=1)

        # find idx nearest corner
        idx_r, idx_c = find_value_rowwise(distances, np.atleast_2d(values).T)
        # spatial relations and nearest corners
        # return array contains only one nearest corner per rect (n, 2=parameter)
        distances = distances[idx_r, idx_c] - self.dot_radii[idx_r]
        directions = np.arctan2(
            xy_dist[idx_r, 1, idx_c],
            xy_dist[idx_r, 0, idx_c])
        return np.array([distances, directions]).T

    def edge_cardinal_relations(self, nearest_edges: bool = True) -> NDArray:
        """Cardinal spatial relations between nearest or farthest edges and dot center

        returns NaN if no edge is cardinal relation to dot center
        """

        edge_points = np.array([(0, 1),  # north
                                (1, 2),  # east
                                (3, 2),  # south
                                (0, 3)  # west
                                ])
        # find nearest edge cross points to dot centers
        # (ecp -> project dot center to rect edge)
        # ecp_xy_diff (n, 2=xy, 4=edges)
        ecp = np.empty_like(self._corners)
        for i in range(4):
            # project of dot center to edge
            ecp[:, :, i] = geometry.line_point_othogonal(
                p1_line=self._corners[:, :, edge_points[i, 0]],
                p2_line=self._corners[:, :, edge_points[i, 1]],
                p3=self.dot_xy,
                outside_segment_nan=True)

        # Euclidean distance of each edge with dot center (n, 4=edges)
        ecp_xy_diff = ecp - np.atleast_3d(self.dot_xy)  # type: ignore
        distances = np.hypot(ecp_xy_diff[:, 0, :],
                             ecp_xy_diff[:, 1, :])
        # find nearest or farthest edge (n, 2=[rect_idx, edge_idx])
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', 'All-NaN slice encountered')
            if nearest_edges:
                values = np.nanmin(distances, axis=1)
            else:
                values = np.nanmax(distances, axis=1)
        # idx_nearest edge
        idx_r, idx_c = find_value_rowwise(distances, np.atleast_2d(values).T)
        # return array contains nearest edge per rect (n, 2=parameter)
        edge_rel = np.full(self.rect_xy.shape, np.nan)
        # distances
        edge_rel[idx_r, 0] = distances[idx_r, idx_c] - self.dot_radii[idx_r]
        # directions
        edge_rel[idx_r, 1] = np.arctan2(
            ecp_xy_diff[idx_r, 1, idx_c],
            ecp_xy_diff[idx_r, 0, idx_c])
        return edge_rel

    def spatial_relations(self, _min: bool = False) -> NDArray:
        """spatial relation between rectangles and dots.

        returns distance long the line of the spatial relation
                and angle  (cached result)
        """
        if self._spatial_relations is None:
            # find closer cardinal edge relation
            corner_rel = self.corner_relations(nearest_corners=_min)
            edge_rel = self.edge_cardinal_relations()
            idx_ccer = corner_rel[:, 0] > edge_rel[:, 0]

            # spat rel shape =(n, 2=parameter)
            self._spatial_relations = np.empty(self.rect_xy.shape)
            self._spatial_relations[idx_ccer, :] = edge_rel[idx_ccer, :]
            self._spatial_relations[~idx_ccer, :] = corner_rel[~idx_ccer, :]

        return self._spatial_relations

    def overlaps(self, minimum_distance: float = 0) -> Union[np.bool_, NDArray[np.bool_]]:
        """Check if rectangles overlap with dots"""
        comp = self.spatial_relations()[:, 0] < minimum_distance
        if comp.shape[0] > 1:
            return comp
        else:
            return comp[0]

    def required_displacements(self, minimum_distance: float = 0) -> NDArray:
        """Minimum required displacement of rectangle B to have no overlap with
        rectangle A.

        Two objects do not overlap if replacements[:, 0] == 0.

        Returns:
            replacements as array of vectors in cartesian coordinates
            if distance (replacements[:, 0]) is 0, no replacement required.

            Calculated movement direction (replacements[:, 0]) is always along the
            line between the two object center.
        """

        sr = self.spatial_relations()
        sr[:, 0] = sr[:, 0] - minimum_distance

        i = np.flatnonzero(sr[:, 0] < 0)  # overlaps
        sr[i, 1] = np.pi + sr[i, 1]  # opposite direction
        # set all nan and override overlapping relations
        rtn = np.full(sr.shape, np.nan)
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
