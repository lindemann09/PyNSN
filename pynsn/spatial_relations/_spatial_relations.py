"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import numpy as np
from numpy.typing import NDArray

from .._lib import geometry

from .._object_arrays.dot_array import BaseDotArray
from .._object_arrays.rectangle_array import BaseRectangleArray
from ._abc_spatial_relations import ABCSpatialRelations


# all init-functions require 2D arrays, all properties are cached
class DotDot(ABCSpatialRelations):

    def __init__(self,
                 dots_a: BaseDotArray,
                 dots_b: BaseDotArray,
                 a_relative_to_b: bool = False) -> None:

        super().__init__(dots_a.xy, dots_b.xy,
                         a_relative_to_b=a_relative_to_b)
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
    def distances_rho(self) -> NDArray:
        if self._distances_rho is None:
            self._distances_rho = np.hypot(self._xy_diff[:, 0],
                                           self._xy_diff[:, 1]) \
                - self._radii_a - self._radii_b
        return self._distances_rho

    @property
    def distances_xy(self) -> NDArray:
        if self._distances_xy is None:
            self._distances_xy = np.abs(self._xy_diff) \
                - np.atleast_2d(self._radii_a + self._radii_b).T
        return self._distances_xy

    def is_inside(self) -> NDArray:
        if self._a_relative_to_b:
            return self.distances_rho < -2 * self._radii_a
        else:
            return self.distances_rho < -2 * self._radii_b

    def displacement_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        return -1*self.distances_rho + minimum_gap


class RectangleRectangle(ABCSpatialRelations):

    def __init__(self,
                 rectangles_a: BaseRectangleArray,
                 rectangles_b: BaseRectangleArray,
                 a_relative_to_b: bool = True) -> None:

        super().__init__(a_xy=rectangles_a.xy, b_xy=rectangles_b.xy,
                         a_relative_to_b=a_relative_to_b)
        self.a_sizes_div2 = rectangles_a.sizes / 2
        self.b_sizes_div2 = rectangles_b.sizes / 2
        a = len(self.a_sizes_div2)
        b = len(self.b_sizes_div2)
        if b == 1:
            self.b_sizes_div2 = self.b_sizes_div2 * np.ones((a, 1))
        elif a == 1:
            self.a_sizes_div2 = self.a_sizes_div2 * np.ones((b, 1))

    @property
    def distances_xy(self) -> NDArray:
        if self._distances_xy is None:
            self._distances_xy = np.abs(self._xy_diff) - \
                (self.a_sizes_div2 + self.b_sizes_div2)

        return self._distances_xy

    def is_inside(self) -> NDArray:
        if self._a_relative_to_b:
            sizes = self.a_sizes_div2
        else:
            sizes = self.b_sizes_div2
        return np.all(self.distances_xy < -1*sizes, axis=1)

    @property
    def distances_rho(self) -> NDArray:

        if self._distances_rho is None:
            # calc distance_inside rect along_axis between object center
            d_a = geometry.center_edge_distance(angles=self.rho,
                                                rect_sizes_div2=self.a_sizes_div2)
            d_b = geometry.center_edge_distance(angles=self.rho,
                                                rect_sizes_div2=self.b_sizes_div2)
            # distance between center minus inside rectangles
            self._distances_rho = np.hypot(self._xy_diff[:, 0],
                                           self._xy_diff[:, 1]) - d_a - d_b

        return self._distances_rho

    def displacement_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        # calc distance_along_line between rect center
        minimum_dist_rho = geometry.center_edge_distance(
            angles=self.rho,
            rect_sizes_div2=np.full(self.a_sizes_div2.shape, minimum_gap))
        return -1*self.distances_rho + minimum_dist_rho


class RectangleDot(ABCSpatialRelations):

    def __init__(self,
                 rectangles: BaseRectangleArray,
                 dots: BaseDotArray,
                 a_relative_to_b: bool = True) -> None:

        self.rect_xy = rectangles.xy
        self.rect_sizes_div2 = rectangles.sizes / 2
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
            self.rect_sizes_div2 = self.rect_sizes_div2 * ones

        assert self.rect_xy.shape == self.dot_xy.shape

        self.__corners = None

        super().__init__(a_xy=self.rect_xy, b_xy=self.dot_xy,
                         a_relative_to_b=a_relative_to_b)

    @property
    def distances_rho(self) -> NDArray:
        if self._distances_rho is None:
            # calc distance_inside rect along_axis between object center
            d_a = geometry.center_edge_distance(angles=self.rho,
                                                rect_sizes_div2=self.rect_sizes_div2)
            # distance between center minus inside rectangles and radii
            self._distances_rho = np.hypot(self._xy_diff[:, 0],
                                           self._xy_diff[:, 1]) \
                - d_a - self.dot_radii

        return self._distances_rho

    @property
    def distances_xy(self) -> NDArray:
        if self._distances_xy is None:
            self._distances_xy = np.abs(self._xy_diff) - self.rect_sizes_div2 \
                - np.atleast_2d(self.dot_radii).T
        return self._distances_xy

    def is_inside(self) -> NDArray:
        if self._a_relative_to_b:
            # rectangle in dots?
            xy_dist_crn = self._corners - \
                np.atleast_3d(self.dot_xy)  # type: ignore
            dot_corner_distances = np.hypot(
                xy_dist_crn[:, 0, :], xy_dist_crn[:, 1, :])  # (n, 4)
            return np.all(dot_corner_distances <= np.atleast_2d(self.dot_radii).T,
                          axis=1)
        else:
            # dots in rectangle
            radii2 = self.dot_radii * np.ones((2, 1))
            # xy_distance + r < rect_size/2
            return np.all(np.abs(self._xy_diff) + radii2.T < self.rect_sizes_div2,
                          axis=1)

    def displacement_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        # Target x and y distance between center to have no overlap at either
        # X or y axes  (TODO this procedure overestimate when close to corner)
        td_xy = self.rect_sizes_div2 + np.atleast_2d(self.dot_radii).T \
            + minimum_gap
        # Target distances along the line between centers, that is, calculate
        # Euclidean distances along the polar radius (rho) that
        # correspond to x and y distances along the cartesian x and y.
        td_center = np.empty_like(td_xy)
        td_center[:, 0] = np.abs(td_xy[:, 0] / np.cos(self.rho))
        td_center[:, 1] = np.abs(td_xy[:, 1] / np.sin(self.rho))
        # find shortest target distance (horizontal or vertical) and subtract
        # the actual distance between center
        return np.min(td_center, axis=1) \
            - np.hypot(self._xy_diff[:, 0], self._xy_diff[:, 1])

    @property
    def _corners(self) -> NDArray:
        """tensor (n, 2, 4) with xy values of the four corners of the rectangles
        0=left-top, 1=right-top, 2=right-bottom, 3=left-bottom
        """
        if self.__corners is None:
            self.__corners = geometry.corners(rect_xy=self.rect_xy,
                                              rect_sizes_div2=self.rect_sizes_div2,
                                              lt_rb_only=False)
        return self.__corners


class DotCoordinate(DotDot):

    def __init__(self,
                 dots: BaseDotArray,
                 coord_xy: NDArray,
                 a_relative_to_b: bool = True) -> None:

        dots_b = BaseDotArray(xy=coord_xy,
                              diameter=np.zeros(coord_xy.shape[0]))
        super().__init__(dots_a=dots, dots_b=dots_b,
                         a_relative_to_b=a_relative_to_b)


class RectangleCoordinate(RectangleRectangle):

    def __init__(self,
                 rectangles: BaseRectangleArray,
                 coord_xy: NDArray,
                 a_relative_to_b: bool = True) -> None:
        rects_b = BaseRectangleArray(xy=coord_xy,
                                     sizes=np.zeros(coord_xy.shape))
        super().__init__(rectangles_a=rectangles, rectangles_b=rects_b,
                         a_relative_to_b=a_relative_to_b)
