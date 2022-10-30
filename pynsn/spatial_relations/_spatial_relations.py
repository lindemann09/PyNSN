"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from cmath import nan
from typing import Union
import numpy as np
from numpy.typing import NDArray

from .._lib import geometry, np_tools
from .._object_arrays import BaseDotArray, BaseRectangleArray
from ._abc_spatial_relations import ABCSpatialRelations


# all init-functions require 2D arrays, all properties are cached
class DotDot(ABCSpatialRelations):

    def __init__(self,
                 a_dots: BaseDotArray,
                 b_dots: BaseDotArray,
                 a_relative_to_b: bool = False) -> None:

        super().__init__(a_array=a_dots,
                         b_array=b_dots,
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
    def distances(self) -> NDArray:
        if self._distances is None:
            self._distances = np.hypot(self._xy_diff[:, 0],
                                       self._xy_diff[:, 1]) \
                - self._a_radii - self._b_radii
        return self._distances

    @property
    def distances_rho(self) -> NDArray:
        return self.distances

    def is_inside(self, minimum_gap: float = 0) -> NDArray:
        if self._a_relative_to_b:
            return self.distances < -2 * self._a_radii - minimum_gap
        else:
            return self.distances < -2 * self._b_radii - minimum_gap

    def _spread_radial_displacement_distances(self, minimum_gap: float) -> NDArray:
        return -1 * self.distances + minimum_gap

    def _spread_displacements(self, minimum_gap: float, polar: bool) -> NDArray:
        displ_polar = np.zeros_like(self._xy_diff)
        i = self.distances < minimum_gap
        displ_polar[i, 0] = -1*self.distances[i] + minimum_gap
        displ_polar[i, 1] = self.rho[i]
        if polar:
            return displ_polar
        else:
            return geometry.polar2cartesian(displ_polar)

    def _gather_displacements(self, minimum_gap: float, polar: bool) -> NDArray:
        """ FIXME not tested"""
        displ_polar = np.zeros_like(self._xy_diff)
        i = self.distances > minimum_gap
        displ_polar[i, 1] = np.pi + self.rho[i]
        if self._a_relative_to_b:
            displ_polar[i, 0] = self.distances + \
                minimum_gap + 2 * self._a_radii[i]
        else:
            displ_polar[i, 1] = self.distances[i] + \
                minimum_gap + 2 * self._b_radii[i]

        if polar:
            return displ_polar
        else:
            return geometry.polar2cartesian(displ_polar)


class RectangleRectangle(ABCSpatialRelations):

    def __init__(self,
                 a_rectangles: BaseRectangleArray,
                 b_rectangles: BaseRectangleArray,
                 a_relative_to_b: bool = True) -> None:

        super().__init__(a_array=a_rectangles,
                         b_array=b_rectangles,
                         a_relative_to_b=a_relative_to_b)
        self.__xy_distances_salted = None
        self.a_sizes_div2 = a_rectangles.sizes / 2
        self.b_sizes_div2 = b_rectangles.sizes / 2
        a = len(self.a_sizes_div2)
        b = len(self.b_sizes_div2)
        if b == 1:
            self.b_sizes_div2 = self.b_sizes_div2 * np.ones((a, 1))
        elif a == 1:
            self.a_sizes_div2 = self.a_sizes_div2 * np.ones((b, 1))

    @property
    def _xy_distances_salted(self):
        """xy_dist salted: xy_distance, randomly salted if x==y ensure one
                preferred cardinal dimension
        """

        if self.__xy_distances_salted is None:
            xy_dist = np.abs(self._xy_diff) - \
                (self.a_sizes_div2 + self.b_sizes_div2)
            # if both dimensions are overlapping identical: choose dimension randomly
            # 'by adding salt'
            i = np.flatnonzero(np.sum(xy_dist < 0, axis=1) == 2)
            xy_dist[i, :] = np_tools.salted_rows(xy_dist[i, :], salt=-1e-30)
            self.__xy_distances_salted = xy_dist

        return self.__xy_distances_salted

    @property
    def distances(self) -> NDArray:
        if self._distances is None:
            xy_dist = self._xy_distances_salted
            n_overlap = np.sum(xy_dist < 0, axis=1)
            self._distances = np.empty(len(xy_dist.shape))

            # 1 dim overlap: max -> distance of not overlapping dimension
            # 2 dim overlap: max -> largest neg. dist. is closest to edge
            i = np.flatnonzero(n_overlap > 0)
            idx = np.argmax(xy_dist[i, :], axis=1)
            self._distances[i] = xy_dist[i, idx]
            # no dimension overlap -> take Euclidean distance between edges
            i = n_overlap == 0
            self._distances[i] = np.hypot(xy_dist[i, 0], xy_dist[i, 1])

        return self._distances

    def _spread_displacements(self, minimum_gap: float, polar: bool) -> NDArray:

        displ_polar = _spread_rectangles(xy_diff=self._xy_diff,
                                         distances=self.distances,
                                         xy_distances_salted=self._xy_distances_salted,
                                         minimum_gap=minimum_gap)
        if self._a_relative_to_b:
            displ_polar[:, 1] = displ_polar[:, 1] + np.pi

        if polar:
            return geometry.polar2cartesian(displ_polar)
        else:
            return displ_polar

    def is_inside(self, minimum_gap: float = 0) -> NDArray:
        if self._a_relative_to_b:
            sizes = self.a_sizes_div2
        else:
            sizes = self.b_sizes_div2
        return np.all(self._xy_distances_salted < -1*sizes - minimum_gap, axis=1)

    def fits_inside(self, minimum_gap: float = 0) -> NDArray:
        if self._a_relative_to_b:
            sml = self.a_sizes_div2 <= self.b_sizes_div2 - minimum_gap
        else:
            sml = self.b_sizes_div2 <= self.a_sizes_div2 - minimum_gap
        return np.all(sml, axis=1)

    @property
    def distances_rho(self) -> NDArray:
        if self._distances_radial is None:
            # calc distance_inside rect along_axis between object center
            d_a = geometry.center_edge_distance(angles=self.rho,
                                                rect_sizes_div2=self.a_sizes_div2)
            d_b = geometry.center_edge_distance(angles=self.rho,
                                                rect_sizes_div2=self.b_sizes_div2)
            # distance between center minus inside rectangles
            self._distances_radial = np.hypot(self._xy_diff[:, 0],
                                              self._xy_diff[:, 1]) - d_a - d_b

        return self._distances_radial

    def _spread_radial_displacement_distances(self, minimum_gap: float) -> NDArray:
        # calc distance_along_line between rect center
        minimum_dist_rho = geometry.center_edge_distance(
            angles=self.rho,
            rect_sizes_div2=np.full(self.a_sizes_div2.shape, minimum_gap))
        return -1*self.distances_rho + minimum_dist_rho

    def _gather_displacements(self, minimum_gap: float, polar: bool) -> NDArray:
        """ FIXME not tested"""

        if self._a_relative_to_b:
            xy_movement = _gather_rectangles(xy_diff=self._xy_diff,
                                             a_sizes_div2=self.b_sizes_div2,
                                             b_sizes_div2=self.a_sizes_div2,
                                             minimum_gap=minimum_gap)
        else:
            xy_movement = _gather_rectangles(xy_diff=self._xy_diff,
                                             a_sizes_div2=self.a_sizes_div2,
                                             b_sizes_div2=self.b_sizes_div2,
                                             minimum_gap=minimum_gap)
        if polar:
            return geometry.cartesian2polar(xy_movement)
        else:
            return xy_movement


class RectangleDot(ABCSpatialRelations):

    def __init__(self,
                 a_rectangles: BaseRectangleArray,
                 b_dots: BaseDotArray,
                 a_relative_to_b: bool = True) -> None:

        super().__init__(a_array=a_rectangles,
                         b_array=b_dots,
                         a_relative_to_b=a_relative_to_b)

        self.__corner_rel = None
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

    @property
    def distances_rho(self) -> NDArray:
        if self._distances_radial is None:
            # calc distance_inside rect along_axis between object center
            d_a = geometry.center_edge_distance(angles=self.rho,
                                                rect_sizes_div2=self.rect_sizes_div2)
            # distance between center minus inside rectangles and radii
            self._distances_radial = np.hypot(self._xy_diff[:, 0],
                                              self._xy_diff[:, 1]) \
                - d_a - self.dot_radii.T  # T-> because it's a column

        return self._distances_radial

    @property
    def distances(self) -> NDArray:
        raise NotImplementedError()

    @property
    def _corner_rel(self) -> NDArray:
        """return distances and directions of all corners to dot edge.
        Array dimensions: (n, 2 [dist, rho], 4 [corners])
        """
        if self.__corner_rel is None:
            # xy dist array of corners : n, 2 [x,y], 4 [corner])
            xy_dist_crn = geometry.corners(rect_xy=self.rect_xy,
                                           rect_sizes_div2=self.rect_sizes_div2,
                                           lt_rb_only=False) \
                - np.atleast_3d(self.dot_xy)  # type: ignore

            dist_corners = np.hypot(xy_dist_crn[:, 0, :], xy_dist_crn[:, 1, :]) \
                - self.dot_radii  # (n_corner, 4)

            # n, 2 [dist, rho], 4 [corners]
            self.__corner_rel = np.full_like(
                xy_dist_crn, fill_value=np.nan)  # all nan
            # distance for corners_outside
            self.__corner_rel[:, 0, :] = dist_corners
            # rho for corners_outside
            self.__corner_rel[:, 1, :] = np.arctan2(xy_dist_crn[:, 1, :],
                                                    xy_dist_crn[:, 0, :])

        return self.__corner_rel

    def is_inside(self, minimum_gap: float = 0) -> NDArray:
        if self._a_relative_to_b:
            # rectangles in dots: distances <= minimum_gap
            return np.all(self._corner_rel[:, 0, :] <= -1*minimum_gap, axis=1)
        else:
            # dots in rectangle
            radii2 = self.dot_radii * np.ones((1, 2))
            # xy_difference + r < rect_size/2 - minimum_gap
            return np.all(np.abs(self._xy_diff) + radii2 <
                          self.rect_sizes_div2 - minimum_gap, axis=1)

    def fits_inside(self, minimum_gap: float = 0) -> NDArray:
        if self._a_relative_to_b:
            # rectangles in dots
            # diagonal of rect fits in circle
            diag = np.atleast_2d(np.hypot(self.rect_sizes_div2[:, 0],
                                          self.rect_sizes_div2[:, 1])).T
            sml = diag <= self.dot_radii - minimum_gap
        else:
            # dots in rectangle
            sml = self.dot_radii * np.ones((1, 2)) \
                <= self.rect_sizes_div2 - minimum_gap
        return np.all(sml, axis=1)

    def _spread_radial_displacement_distances(self, minimum_gap: float) -> NDArray:
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

    def _spread_displacements(self, minimum_gap: float, polar: bool) -> NDArray:
        if self._a_relative_to_b:
            # corner relations to polar (n_corner, xy)
            corner_rel = self._corner_rel
            # find closest corner
            min_dist = np.atleast_2d(np.nanmin(corner_rel[:, 0, :], axis=1)).T
            i_r, i_c = np_tools.find_rowwise_one_true(
                corner_rel[:, 0, :] == min_dist)

            displ_polar = np.zeros_like(self._xy_diff)
            displ_polar[i_r, :] = corner_rel[i_r, :, i_c]
            # adapt distances
            displ_polar[:, 0] = -1*displ_polar[:, 0] + minimum_gap
        else:
            radii2 = self.dot_radii * np.ones((1, 2))
            xy_dist = np.abs(self._xy_diff) - (self.rect_sizes_div2 + radii2)
            # if both dimensions are overlapping identical: choose dimension randomly
            # 'by adding salt'
            i = np.flatnonzero(np.sum(xy_dist < 0, axis=1) == 2)
            xy_dist[i, :] = np_tools.salted_rows(xy_dist[i, :], salt=-1e-30)
            displ_polar = _spread_rectangles(xy_diff=self._xy_diff,
                                             distances=self.distances,
                                             xy_distances_salted=xy_dist,
                                             minimum_gap=minimum_gap)

        if polar:
            return geometry.polar2cartesian(displ_polar)
        else:
            return displ_polar

    def _gather_displacements(self, minimum_gap: float, polar: bool) -> NDArray:
        """ """
        if self._a_relative_to_b:
            # corner relations to polar (n_corner, xy)
            corner_rel = self._corner_rel
            # find those requiring a replacement
            i_r, i_c = np.nonzero(corner_rel[:, 0, :] > -1*minimum_gap)
            # polar-to-cartesian for each required corner movement
            xy_move_corner = np.zeros_like(corner_rel)
            xy_move_corner_all = np.zeros_like(corner_rel)
            # move always in opposite direct
            diag = np.atleast_2d(np.hypot(self.rect_sizes_div2[:, 0],
                                          self.rect_sizes_div2[:, 1])).T
            rho = np.pi + corner_rel[i_r, 1, i_c]
            dist = corner_rel[i_r, 0, i_c] + minimum_gap

            xy_move_corner[i_r, 0, i_c] = dist * np.cos(rho)
            xy_move_corner[i_r, 1, i_c] = dist * np.sin(rho)
            # find largest x and y movement required for each rect
            #   -> abs minimum across corners: (n, 2):
            print(xy_move_corner[:, 1, 0:4])
            xy_movement = np_tools.abs_maximum(xy_move_corner, axis=2)
        else:
            radii2 = self.dot_radii * np.ones((1, 2))
            xy_movement = _gather_rectangles(xy_diff=self._xy_diff,
                                             a_sizes_div2=self.rect_sizes_div2,
                                             b_sizes_div2=radii2,
                                             minimum_gap=minimum_gap)
        if polar:
            return geometry.cartesian2polar(xy_movement)
        else:
            return xy_movement


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


def _spread_rectangles(xy_diff: NDArray,
                       distances: NDArray,
                       xy_distances_salted: NDArray,
                       minimum_gap: float) -> NDArray:
    # 1 dim overlap: max -> distance of not overlapping dimension
    # 2 dim overlap: max -> largest neg. dist. is closest to edge
    cardinal_relations = np.array([
        np.where(xy_diff[:, 0] > 0, 0, np.pi),
        np.where(xy_diff[:, 1] > 0,
                 geometry.NORTH, geometry.SOUTH)
    ]).T

    i = np.flatnonzero(distances < minimum_gap)
    idx = np.argmax(xy_distances_salted[i, :], axis=1)

    displ_polar = np.zeros_like((len(xy_diff), 2))
    displ_polar[i, 0] = -1*distances[i] + minimum_gap
    displ_polar[i, 1] = cardinal_relations[i, idx]
    return displ_polar


def _gather_rectangles(xy_diff: NDArray,
                       a_sizes_div2:  NDArray,
                       b_sizes_div2: NDArray,
                       minimum_gap: float) -> NDArray:
    """helper: calculates the required displacement to move rect b into rect a
    """
    # xy_movement: neg numbers indicate no movement required, that is,
    #   there already full overlap on that dimension
    xy_movement = np.abs(xy_diff) - (
        a_sizes_div2 - b_sizes_div2 - minimum_gap)  # type: ignore
    # remove negative distances
    xy_movement[np.where(xy_movement <= 0)] = 0
    # adjust direction
    return xy_movement * np.sign(xy_diff) * -1


def relations(a_array: Union[NDArray, BaseDotArray, BaseRectangleArray],
              b_array: Union[NDArray, BaseDotArray, BaseRectangleArray],
              a_relative_to_b: bool = False) -> ABCSpatialRelations:
    """returns the required objects to calculate the spatial relations between
    array_a and array_b.
    """
    if isinstance(a_array, BaseDotArray):
        if isinstance(b_array, BaseDotArray):
            return DotDot(a_array, b_array, a_relative_to_b)
        elif isinstance(b_array, BaseRectangleArray):
            return RectangleDot(b_array, a_array, not a_relative_to_b)
        elif isinstance(b_array, np.ndarray):
            return DotCoordinate(a_array, b_array, a_relative_to_b)
    elif isinstance(a_array, BaseRectangleArray):
        if isinstance(b_array, BaseDotArray):
            return RectangleDot(a_array, b_array, a_relative_to_b)
        elif isinstance(b_array, BaseRectangleArray):
            return RectangleRectangle(b_array, a_array, a_relative_to_b)
        elif isinstance(b_array, np.ndarray):
            return RectangleCoordinate(a_array, b_array, a_relative_to_b)
    elif isinstance(a_array, np.ndarray):
        if isinstance(b_array, BaseDotArray):
            return DotCoordinate(b_array, a_array, not a_relative_to_b)
        elif isinstance(b_array, BaseRectangleArray):
            return RectangleCoordinate(b_array, a_array, not a_relative_to_b)

    raise TypeError(f"Spatial relations class for {type(a_array)} and "
                    f"{type(a_array)} does not exist.")
