"""Numpy optimized  geometry functions for the relations of Dots and Rectangles"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

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
        if self._cache_distances is None:
            self._cache_distances = np.hypot(self._xy_diff[:, 0],
                                             self._xy_diff[:, 1]) \
                - self._a_radii - self._b_radii
        return self._cache_distances

    @property
    def distances_radial(self) -> NDArray:
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
            displ_polar[i, 0] = self.distances[i] + \
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
        self._cache_xy_distances_salted = None
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

        if self._cache_xy_distances_salted is None:
            xy_dist = np.abs(self._xy_diff) - \
                (self.a_sizes_div2 + self.b_sizes_div2)
            # if both dimensions are overlapping identical: choose dimension randomly
            # 'by adding salt'
            i = np.flatnonzero(np.sum(xy_dist < 0, axis=1) == 2)
            xy_dist[i, :] = np_tools.salted_rows(xy_dist[i, :], salt=-1e-30)
            self._cache_xy_distances_salted = xy_dist

        return self._cache_xy_distances_salted

    @property
    def distances(self) -> NDArray:
        if self._cache_distances is None:
            xy_dist = self._xy_distances_salted
            n_overlap = np.sum(xy_dist < 0, axis=1)
            self._cache_distances = np.empty(len(xy_dist))

            # 1 dim overlap: max -> distance of not overlapping dimension
            # 2 dim overlap: max -> largest neg. dist. is closest to edge
            i = np.flatnonzero(n_overlap > 0)
            idx = np.argmax(xy_dist[i, :], axis=1)
            self._cache_distances[i] = xy_dist[i, idx]
            # no dimension overlap -> take Euclidean distance between edges
            i = n_overlap == 0
            self._cache_distances[i] = np.hypot(
                xy_dist[i, 0], xy_dist[i, 1])

        return self._cache_distances

    @property
    def distances_radial(self) -> NDArray:
        if self._cache_distances_radial is None:
            # calc distance_inside rect along_axis between object center
            d_a = geometry.center_edge_distance(angles=self.rho,
                                                rect_sizes_div2=self.a_sizes_div2)
            d_b = geometry.center_edge_distance(angles=self.rho,
                                                rect_sizes_div2=self.b_sizes_div2)
            # distance between center minus inside rectangles
            self._cache_distances_radial = np.hypot(self._xy_diff[:, 0],
                                                    self._xy_diff[:, 1]) - d_a - d_b

        return self._cache_distances_radial

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

    def _spread_radial_displacement_distances(self, minimum_gap: float) -> NDArray:
        # calc distance_along_line between rect center
        minimum_dist_rho = geometry.center_edge_distance(
            angles=self.rho,
            rect_sizes_div2=np.full(self.a_sizes_div2.shape, minimum_gap))
        return -1*self.distances_radial + minimum_dist_rho

    def _spread_displacements(self, minimum_gap: float, polar: bool) -> NDArray:

        displ_polar = _spread_rectangles(xy_diff=self._xy_diff,
                                         distances=self.distances,
                                         xy_dist=self._xy_distances_salted,
                                         minimum_gap=minimum_gap)
        if self._a_relative_to_b:
            displ_polar[:, 1] = displ_polar[:, 1] + np.pi

        if polar:
            return displ_polar
        else:
            return geometry.polar2cartesian(displ_polar)

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

        self._cache_corner_rel = None
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
    def distances_radial(self) -> NDArray:
        if self._cache_distances_radial is None:
            # calc distance_inside rect along_axis between object center
            d_a = geometry.center_edge_distance(angles=self.rho,
                                                rect_sizes_div2=self.rect_sizes_div2)
            # distance between center minus inside rectangles and radii
            self._cache_distances_radial = np.hypot(self._xy_diff[:, 0],
                                                    self._xy_diff[:, 1]) \
                - d_a - self.dot_radii.T  # T-> because it's a column

        return self._cache_distances_radial

    @property
    def distances(self) -> NDArray:
        if self._cache_distances is None:
            corner_rel = self._corner_relations(angles=False)
            radii2 = self.dot_radii * np.ones((1, 2))
            # xy_dist: if not overlapping -> inf.
            xy_dist = np.full_like(self._xy_diff, fill_value=np.inf)
            # rectangles that overlap with center of the dot
            overlaps = np.abs(self._xy_diff) < self.rect_sizes_div2  # overlaps
            # get the dimensions of overlapping one
            idx = np.nonzero(overlaps[:, (1, 0)])
            xy_dist[idx] = np.abs(self._xy_diff[idx]) \
                - self.rect_sizes_div2[idx] - radii2[idx]

            # join Euclidean corner distance and cardinal distances of overlaps
            all_dist = np.hstack([corner_rel[:, 0, :], xy_dist])

            self._cache_distances = np.min(all_dist, axis=1)

        return self._cache_distances

    def _corner_relations(self, angles=True) -> NDArray:
        """return radial distances and directions(angle) of all corners to dot edge.
        Array dimensions: (n, 2 [dist, angle], 4 [corners])

        negative distance indicate corner inside dot

        if angles=False, angles will not be calculated (for efficiently)
        """
        # angles desired and required,
        # because angles missing or nothing calculated so far?
        calc_angles = angles and \
            ((self._cache_corner_rel is not None and
              np.isnan(self._cache_corner_rel[0, 1, 0])) or
             self._cache_corner_rel is None)

        if self._cache_corner_rel is None or calc_angles:
            # xy dist array of corners : n, 2 [x,y], 4 [corner])
            xy_dist_crn = geometry.corners(rect_xy=self.rect_xy,
                                           rect_sizes_div2=self.rect_sizes_div2,
                                           lt_rb_only=False) \
                - np.atleast_3d(self.dot_xy)  # type: ignore

            if self._cache_corner_rel is None:
                # make array : n, 2 [dist, rho], 4 [corners]
                self._cache_corner_rel = np.full_like(xy_dist_crn,
                                                      fill_value=np.nan)
                # distance for corners
                self._cache_corner_rel[:, 0, :] = np.hypot(
                    xy_dist_crn[:, 0, :], xy_dist_crn[:, 1, :])\
                    - self.dot_radii  # (n_corner, 4)

            if calc_angles:
                # angle for corners_outside
                self._cache_corner_rel[:, 1, :] = np.arctan2(xy_dist_crn[:, 1, :],
                                                             xy_dist_crn[:, 0, :])

        return self._cache_corner_rel

    def is_inside(self, minimum_gap: float = 0) -> NDArray:
        if self._a_relative_to_b:
            # rectangles in dots: distances <= minimum_gap
            corner_rel = self._corner_relations(angles=False)
            return np.all(corner_rel[:, 0, :] <= -1*minimum_gap, axis=1)
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
            # get corner
            corner = geometry.corners(rect_xy=self.rect_xy,
                                      rect_sizes_div2=self.rect_sizes_div2,
                                      lt_rb_only=False)
            # get relations and set corners outside to nan
            idx_rect, idx_corner = np.nonzero(
                self._corner_relations(angles=False)[:, 0, :] > minimum_gap)
            corner[idx_rect, :, idx_corner] = np.nan

            intersections = [geometry.line_circle_intersection(
                line_points=corner[:, :, x],
                line_directions=self.rho,
                circle_center=self.dot_xy,
                circle_radii=self.dot_radii) for x in range(4)]
            intersections = np.stack(intersections, axis=2)

            # get idx of corner with maximum distance
            dist_xy = intersections - corner
            distances = np.hypot(dist_xy[:, 0, :],
                                 dist_xy[:, 1, :])  # distance to intersection

            idx_rect, idx_corner = np_tools.nonzero_one_per_dim(
                distances == np.atleast_2d(np.nanmax(distances, axis=1)).T)

            # get respective intersections and calc displacement (cartesian)
            displ_cart = np.zeros_like(self.rect_sizes_div2)
            displ_cart[idx_rect, :] = intersections[idx_rect, :, idx_corner] - \
                corner[idx_rect, :, idx_corner]

            if polar:
                return geometry.cartesian2polar(displ_cart)
            else:
                return displ_cart

        else:
            radii2 = self.dot_radii * np.ones((1, 2))
            xy_dist = np.abs(self._xy_diff) - (self.rect_sizes_div2 + radii2)
            displ_polar = _spread_rectangles(xy_diff=self._xy_diff,
                                             distances=self.distances,
                                             xy_dist=xy_dist,
                                             minimum_gap=minimum_gap)
            if polar:
                return displ_polar
            else:
                return geometry.polar2cartesian(displ_polar)

    def _gather_displacements(self, minimum_gap: float, polar: bool) -> NDArray:
        """ """
        if self._a_relative_to_b:
            # distance points circle edge

            corner_rel = self._corner_relations(angles=True)
            # find those requiring a replacement
            i_r, i_c = np.nonzero(corner_rel[:, 0, :] > -1*minimum_gap)
            # polar-to-cartesian for each required corner movement
            rho = np.pi + corner_rel[i_r, 1, i_c]
            dist = minimum_gap + corner_rel[i_r, 0, i_c]
            xy_move_corner = np.zeros_like(corner_rel)
            xy_move_corner[i_r, 0, i_c] = dist * np.cos(rho)
            xy_move_corner[i_r, 1, i_c] = dist * np.sin(rho)

            # find largest x and y movement required for each rect
            #   -> abs minimum across corners: (n, 2):
            xy_movement = np_tools.abs_maximum(xy_move_corner, axis=2)
            # Find rectangles that are still not fully inside circle and
            # place them to the out possible position at the left, right, top or
            # bottom of the circle.
            # Explanation: Some are placed to far from center and stick out..
            # This are the dots that overlap with dot center at one dimension.
            # Solution: Find max_xy_diff (center differences) by calculating
            # required chord position and replace these rectangles
            #
            # swap xy column for chord len, because width calculates y and height calculates x distance
            max_xy_diff = geometry.chord_distance_to_center(
                r=self.dot_radii,
                chord_len=2*self.rect_sizes_div2[:, (1, 0)]) - self.rect_sizes_div2  # type: ignore
            # find coordinates those too far, by comparing max_xy_diff
            # with abs of new position
            # take candidates that have at least one  too far (i),
            i = np.any(np.abs(self._xy_diff + xy_movement) > max_xy_diff,
                       axis=1)
            # target xy_diff is maximum diff,
            target_xy_diff = max_xy_diff[i, :]
            # take target_xy_diff is corresponds to large xy_movement,
            # that is, set min(move_set) to zero (take randomly one in case of
            # identical values)
            move_size = np.abs(xy_movement[i, :])
            idx = np_tools.nonzero_one_per_dim(
                move_size == np.atleast_2d(np.min(move_size, axis=1)).T)
            target_xy_diff[idx] = 0
            # keep side by keeping previous sign,
            xy_movement[i] = np_tools.sign(self._xy_diff[i]) * target_xy_diff \
                - self._xy_diff[i]
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
                       xy_dist: NDArray,
                       minimum_gap: float) -> NDArray:
    """helper: move out in one cardinal axes"""

    rtn = np.zeros_like(xy_diff)
    # if both dimensions are overlapping to an identical extend:
    # choose dimension randomly
    # 'by adding salt'
    i = (xy_dist[:, 0] < 0) & (xy_dist[:, 0] == xy_dist[:, 1])
    xy_dist[i, :] = np_tools.salted_rows(xy_dist[i, :], salt=-1e-30)

    i = np.flatnonzero(distances < minimum_gap)
    cardinal_relations = np.array([
        np.where(xy_diff[i, 0] < 0, 0, np.pi),
        np.where(xy_diff[i, 1] < 0, geometry.NORTH, geometry.SOUTH)]).T
    all_col = range(cardinal_relations.shape[0])  # all columns
    # where is the max displacement
    which_coord = np.argmax(xy_dist[i, :], axis=1)
    rtn[i, 0] = minimum_gap - xy_dist[i, which_coord]
    rtn[i, 1] = cardinal_relations[all_col, which_coord]
    return rtn


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
    return xy_movement * np_tools.sign(xy_diff) * -1


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
