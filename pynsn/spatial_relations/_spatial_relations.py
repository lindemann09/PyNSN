"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Union

import warnings
import numpy as np
from numpy.typing import NDArray

from .._lib.np_tools import find_value_rowwise
from .._lib import geometry

from .._object_arrays.dot_array import BaseDotArray
from .._object_arrays.rectangle_array import BaseRectangleArray
from ._abc_spatial_relations import ABCSpatialRelations


# all init-functions require 2D arrays
class CoordinateCoordinate(ABCSpatialRelations):

    def __init__(self,
                 a_xy: NDArray[np.floating],
                 b_xy: NDArray[np.floating],
                 A_relative_to_B: bool = False) -> None:
        super().__init__(a_xy, b_xy, A_relative_to_B=A_relative_to_B)

    def is_inside(self) -> NDArray:
        return np.full(len(self._xy_diff), np.nan)

    @ property
    def distances(self) -> NDArray:
        if self._distances is None:
            self._distances = np.hypot(
                self._xy_diff[:, 0], self._xy_diff[:, 1])
        return self._distances

    def displacement_distances(self, minimum_distance: float = 0) -> NDArray:
        raise NotImplementedError()


class DotDot(ABCSpatialRelations):

    def __init__(self,
                 dots_a: BaseDotArray,
                 dots_b: BaseDotArray,
                 A_relative_to_B: bool = False) -> None:

        super().__init__(dots_a.xy, dots_b.xy,
                         A_relative_to_B=A_relative_to_B)
        self._radii_a = dots_a.diameter / 2
        self._radii_b = dots_b.diameter / 2
        a = len(self._radii_a)
        b = len(self._radii_b)
        if b == 1:
            self._radii_b = self._radii_b * np.ones(a)
        elif a == 1:
            self._radii_a = self._radii_a * np.ones(b)

        assert self._radii_a.shape == self._radii_b.shape

    @property
    def distances(self) -> NDArray:
        if self._distances is None:
            self._distances = np.hypot(self._xy_diff[:, 0],
                                       self._xy_diff[:, 1]) \
                - self._radii_a - self._radii_b
        return self._distances

    def is_inside(self) -> NDArray:
        if self._A_relative_to_B:
            return self.distances < -2 * self._radii_a
        else:
            return self.distances < -2 * self._radii_b

    def displacement_distances(self, minimum_distance: float = 0) -> NDArray:
        # FIXME not tested
        return self.distances - minimum_distance


class RectangleRectangle(ABCSpatialRelations):

    def __init__(self,
                 rectangles_a: BaseRectangleArray,
                 rectangles_b: BaseRectangleArray,
                 A_relative_to_B: bool = True) -> None:

        super().__init__(a_xy=rectangles_a.xy, b_xy=rectangles_b.xy,
                         A_relative_to_B=A_relative_to_B)
        self.a_sizes = rectangles_a.sizes
        self.b_sizes = rectangles_b.sizes
        a = len(self.a_sizes)
        b = len(self.b_sizes)
        if b == 1:
            self.b_sizes = self.b_sizes * np.ones((a, 1))
        elif a == 1:
            self.a_sizes = self.a_sizes * np.ones((b, 1))
        self._xy_diff_rect = np.abs(self._xy_diff) - \
            (self.a_sizes + self.b_sizes) / 2

    def is_inside(self) -> NDArray:
        if self._A_relative_to_B:
            sizes = self.a_sizes
        else:
            sizes = self.b_sizes
        return np.all(self._xy_diff_rect < -1*sizes, axis=1)

    @property
    def distances(self) -> NDArray:
        if self._distances is None:
            xy_diff = np.copy(self._xy_diff_rect)

            # find rows with both coordinate positive or negative (i_pn or i_bn)
            cnt_neg = np.sum(xy_diff < 0, axis=1)
            i_bn = np.flatnonzero(cnt_neg == 2)

            # make the large (closer to zero) coordinate positive
            lrg = xy_diff[i_bn, 0] >= xy_diff[i_bn, 1]
            xy_diff[i_bn[lrg], 0] = -1 * xy_diff[i_bn[lrg], 0]
            xy_diff[i_bn[~lrg], 1] = -1 * xy_diff[i_bn[~lrg], 1]
            # and set the other coordinate to zero
            xy_diff[i_bn[lrg], 1] = 0
            xy_diff[i_bn[~lrg], 0] = 0
            # set all remaining negative values to zero, because
            # the other coordinate are positive and thus define the distance
            xy_diff[np.where(xy_diff < 0)] = 0
            self._distances = np.hypot(xy_diff[:, 0],
                                       xy_diff[:, 1])
            # set distance with previous both negative to negative
            self._distances[i_bn] = self._distances[i_bn] * -1

        return self._distances

    def displacement_distances(self, minimum_distance: float = 0) -> NDArray:
        # calc distance_along_line between rect center

        # but enlarge one rect by minimum distance
        d_a = geometry.center_edge_distance(angles=self.angle,
                                            rect_sizes=self.a_sizes + 2 * minimum_distance)
        d_b = geometry.center_edge_distance(angles=self.angle,
                                            rect_sizes=self.b_sizes)
        # distance between center minus inside rectangles
        return np.hypot(self._xy_diff[:, 0],
                        self._xy_diff[:, 1]) - d_a - d_b


class RectangleDot(ABCSpatialRelations):

    def __init__(self,
                 rectangles: BaseRectangleArray,
                 dots: BaseDotArray,
                 A_relative_to_B: bool = True) -> None:

        self.rect_xy = rectangles.xy
        self.rect_sizes = rectangles.sizes
        self.dot_xy = dots.xy
        self.dot_radii = dots.diameter / 2

        n_rects = len(self.rect_xy)
        n_dots = len(self.dot_xy)
        if n_dots == 1:
            self.dot_xy = self.dot_xy * np.ones((n_rects, 1))
            self.dot_radii = self.dot_radii * np.ones(n_rects)
        elif n_rects == 1:
            ones = np.ones((n_dots, 1))
            self.rect_xy = self.rect_xy * ones
            self.rect_sizes = self.rect_sizes * ones

        assert self.rect_xy.shape == self.dot_xy.shape

        self.__corners = None

        super().__init__(a_xy=self.rect_xy, b_xy=self.dot_xy,
                         A_relative_to_B=A_relative_to_B)

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

    @property
    def distances(self) -> NDArray:
        return NotImplemented

    def is_inside(self) -> NDArray:
        return NotImplemented

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

    def displacement_distances(self, minimum_distance: float = 0) -> NDArray:
        # spatial relations
        corner_rel = self.corner_relations(nearest_corners=False)
        edge_rel = self.edge_cardinal_relations()
        idx_ccer = corner_rel[:, 0] > edge_rel[:, 0]

        # spat rel shape =(n, 2=parameter)
        spatrel = np.empty(self.rect_xy.shape)
        spatrel[idx_ccer, :] = edge_rel[idx_ccer, :]
        spatrel[~idx_ccer, :] = corner_rel[~idx_ccer, :]

        spatrel[:, 0] = spatrel[:, 0] - minimum_distance

        return spat_rel2displacement(spatrel)


class DotCoordinate(DotDot):

    def __init__(self,
                 dots: BaseDotArray,
                 coord_xy: NDArray,
                 A_relative_to_B: bool = True) -> None:

        dots_b = BaseDotArray(xy=coord_xy,
                              diameter=np.zeros(coord_xy.shape[0]))
        super().__init__(dots_a=dots, dots_b=dots_b,
                         A_relative_to_B=A_relative_to_B)


class RectangleCoordinate(RectangleRectangle):

    def __init__(self,
                 rectangles: BaseRectangleArray,
                 coord_xy: NDArray,
                 A_relative_to_B: bool = True) -> None:
        rects_b = BaseRectangleArray(xy=coord_xy,
                                     sizes=np.zeros(coord_xy.shape))
        super().__init__(rectangles_a=rectangles, rectangles_b=rects_b,
                         A_relative_to_B=A_relative_to_B)


# SpatRelArrayType = Union[NDArray, BaseDotArray, BaseRectangleArray]


# def spatial_relations(array_A: SpatRelArrayType,
#                       array_B: SpatRelArrayType) -> ABCSpatialRelations:
#     if isinstance(array_A, BaseRectangleArray):
#         if isinstance(array_B, BaseRectangleArray):
#             return RectangleRectangle(array_A, array_B)
#         elif isinstance(array_B, BaseDotArray):
#             return RectangleDot(array_A, array_B)
#         elif isinstance(array_B, np.ndarray):
#             return RectangleCoordinate(array_A, array_B)

#     elif isinstance(array_A, BaseDotArray):
#         if isinstance(array_B, BaseDotArray):
#             return DotDot(array_A, array_B)
#         elif isinstance(array_B, BaseRectangleArray):
#             return RectangleDot(array_B, array_A, A_relative_to_B=False)
#         elif isinstance(array_B, np.ndarray):
#             return DotCoordinate(array_A, array_B)

#     elif isinstance(array_A, np.ndarray):
#         if isinstance(array_B, np.ndarray):
#             return CoordinateCoordinate(array_A, array_B)
#         if isinstance(array_B, BaseDotArray):
#             return DotCoordinate(array_B, array_A, A_relative_to_B=False)
#         elif isinstance(array_B, BaseRectangleArray):
#             return RectangleCoordinate(array_B, array_A, A_relative_to_B=False)

#     raise NotImplementedError()
