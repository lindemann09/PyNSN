"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from tkinter import dialog
import numpy as np
from numpy.typing import NDArray

from .._lib import geometry

from .._object_arrays.dot_array import BaseDotArray
from .._object_arrays.rectangle_array import BaseRectangleArray
from ._abc_spatial_relations import ABCSpatialRelations


# all init-functions require 2D arrays, all properties are cached
class DotDot(ABCSpatialRelations):

    def __init__(self,
                 a_dots: BaseDotArray,
                 b_dots: BaseDotArray,
                 a_relative_to_b: bool = False) -> None:

        super().__init__(a_dots.xy, b_dots.xy,
                         a_relative_to_b=a_relative_to_b)
        self._a_radii = a_dots.diameter / 2
        self._b_radii = b_dots.diameter / 2
        a = len(self._a_radii)
        b = len(self._b_radii)
        if b == 1:
            self._b_radii = self._b_radii * np.ones(a)
        elif a == 1:
            self._a_radii = self._a_radii * np.ones(b)

        assert self._a_radii.shape == self._b_radii.shape

    def fits_inside(self, minimum_gap: float = 0) -> NDArray:
        if self._a_relative_to_b:
            return self._a_radii <= self._b_radii - minimum_gap
        else:
            return self._b_radii <= self._a_radii - minimum_gap

    @property
    def distances_rho(self) -> NDArray:
        if self._distances_rho is None:
            self._distances_rho = np.hypot(self._xy_diff[:, 0],
                                           self._xy_diff[:, 1]) \
                - self._a_radii - self._b_radii
        return self._distances_rho

    @property
    def distances_xy(self) -> NDArray:
        if self._distances_xy is None:
            self._distances_xy = np.abs(self._xy_diff) \
                - np.atleast_2d(self._a_radii + self._b_radii).T
        return self._distances_xy

    def is_inside(self) -> NDArray:
        if self._a_relative_to_b:
            return self.distances_rho < -2 * self._a_radii
        else:
            return self.distances_rho < -2 * self._b_radii

    def _spread_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        return -1*self.distances_rho + minimum_gap


class RectangleRectangle(ABCSpatialRelations):

    def __init__(self,
                 a_rectangles: BaseRectangleArray,
                 b_rectangles: BaseRectangleArray,
                 a_relative_to_b: bool = True) -> None:

        super().__init__(a_xy=a_rectangles.xy, b_xy=b_rectangles.xy,
                         a_relative_to_b=a_relative_to_b)
        self.a_sizes_div2 = a_rectangles.sizes / 2
        self.b_sizes_div2 = b_rectangles.sizes / 2
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

    def fits_inside(self, minimum_gap: float = 0) -> NDArray:
        if self._a_relative_to_b:
            sml = self.a_sizes_div2 <= self.b_sizes_div2 - minimum_gap
        else:
            sml = self.b_sizes_div2 <= self.a_sizes_div2 - minimum_gap
        return np.all(sml, axis=1)

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

    def _spread_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        # calc distance_along_line between rect center
        minimum_dist_rho = geometry.center_edge_distance(
            angles=self.rho,
            rect_sizes_div2=np.full(self.a_sizes_div2.shape, minimum_gap))
        return -1*self.distances_rho + minimum_dist_rho


class RectangleDot(ABCSpatialRelations):

    def __init__(self,
                 a_rectangles: BaseRectangleArray,
                 b_dots: BaseDotArray,
                 a_relative_to_b: bool = True) -> None:

        self.rect_xy = a_rectangles.xy
        self.rect_sizes_div2 = a_rectangles.sizes / 2
        self.dot_xy = b_dots.xy
        dot_radii = b_dots.diameter / 2

        n_rects = len(self.rect_xy)
        n_dots = len(self.dot_xy)
        if n_dots == 1:
            self.dot_xy = self.dot_xy * np.ones((n_rects, 1))
            dot_radii = dot_radii * np.ones(n_rects)
        elif n_rects == 1:
            ones = np.ones((n_dots, 1))
            self.rect_xy = self.rect_xy * ones
            self.rect_sizes_div2 = self.rect_sizes_div2 * ones

        self.dot_radii = np.atleast_2d(dot_radii).T  # dot radii_column

        assert self.rect_xy.shape == self.dot_xy.shape

        super().__init__(a_xy=self.rect_xy,
                         b_xy=self.dot_xy,
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
                - d_a - self.dot_radii.T  # T-> because it's a column

        return self._distances_rho

    @property
    def distances_xy(self) -> NDArray:
        if self._distances_xy is None:
            self._distances_xy = np.abs(self._xy_diff) - self.rect_sizes_div2 \
                - self.dot_radii
        return self._distances_xy

    def is_inside(self) -> NDArray:
        if self._a_relative_to_b:
            # rectangles in dots
            xy_dist_crn = geometry.corners(rect_xy=self.rect_xy,
                                           rect_sizes_div2=self.rect_sizes_div2,
                                           lt_rb_only=False) \
                - np.atleast_3d(self.dot_xy)  # type: ignore
            dot_corner_distances = np.hypot(xy_dist_crn[:, 0, :],
                                            xy_dist_crn[:, 1, :])  # (n, 4)
            return np.all(dot_corner_distances <= self.dot_radii, axis=1)

        else:
            # dots in rectangle
            radii2 = self.dot_radii * np.ones((1, 2))
            # xy_distance + r < rect_size/2
            return np.all(np.abs(self._xy_diff) + radii2.T < self.rect_sizes_div2,
                          axis=1)

    def fits_inside(self, minimum_gap: float = 0) -> NDArray:
        if self._a_relative_to_b:
            # rectangles in dots
            # diagonal of rect fits in circle
            sml = np.hypot(self.rect_sizes_div2[:, 0], self.rect_sizes_div2[:, 1]) \
                <= self.dot_radii - minimum_gap
        else:
            # dots in rectangle
            sml = self.dot_radii * np.ones((1, 2)) \
                <= self.rect_sizes_div2 - minimum_gap
        return np.all(sml, axis=1)

    def _spread_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        # Target x and y distance between center to have no overlap at either
        # X or y axes  (TODO this procedure overestimate when close to corner)
        td_xy = self.rect_sizes_div2 + self.dot_radii + minimum_gap
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

    def _gather_distances_rho(self, minimum_gap: float = 0) -> NDArray:
        return super()._gather_distances_rho(minimum_gap)

    def _gather_polar_shortest(self, minimum_gap: float = 0) -> NDArray:
        raise NotImplementedError()


class DotCoordinate(DotDot):

    def __init__(self,
                 a_dots: BaseDotArray,
                 b_coord_xy: NDArray,
                 a_relative_to_b: bool = True) -> None:

        dots_b = BaseDotArray(xy=b_coord_xy,
                              diameter=np.zeros(b_coord_xy.shape[0]))
        super().__init__(a_dots=a_dots, b_dots=dots_b,
                         a_relative_to_b=a_relative_to_b)


class RectangleCoordinate(RectangleRectangle):

    def __init__(self,
                 a_rectangles: BaseRectangleArray,
                 b_coord_xy: NDArray,
                 a_relative_to_b: bool = True) -> None:
        rects_b = BaseRectangleArray(xy=b_coord_xy,
                                     sizes=np.zeros(b_coord_xy.shape))
        super().__init__(a_rectangles=a_rectangles, b_rectangles=rects_b,
                         a_relative_to_b=a_relative_to_b)
