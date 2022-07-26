from __future__ import annotations
from typing import Iterable, Iterator

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from abc import ABCMeta, abstractmethod
from copy import deepcopy

import numpy as np
from scipy import spatial

from .point_array import PointArray
from .tools import BrownianMotion
from .._lib import constants
from .._lib import geometry
from .._lib import misc
from .._lib import rng
from .._lib.coordinate import Coordinate
from .._lib.exception import NoSolutionError
from .._lib.lib_typing import List, Union, Tuple, Sequence, NDArray
from .._shapes.dot import Dot
from .._shapes.point import Point
from .._shapes.rectangle import Rectangle


class ABCObjectArray(PointArray, metaclass=ABCMeta):

    @property
    @abstractmethod
    def surface_areas(self):
        """Vector with the surface areas of all objects"""
        pass

    @property
    @abstractmethod
    def perimeter(self) -> NDArray:
        """Vector with the perimeter of all objects"""
        pass

    @abstractmethod
    def as_dict(self):
        """Dict representation of the array"""
        return super().as_dict()

    @staticmethod
    @abstractmethod
    def read_from_dict(the_dict):
        pass

    @staticmethod
    @abstractmethod
    def load(json_file_name):
        pass

    @abstractmethod
    def copy(self, indices: Union[int, Sequence[int], None] =None,
             deep_copy=True) -> ABCObjectArray:
        pass

    @abstractmethod
    def iter_objects(self, indices=None) -> Iterator:
        pass

    @abstractmethod
    def add(self, obj):
        return super().add(obj)

    @abstractmethod
    def find_objects(self, size=None, attribute=None):
        pass

    @abstractmethod
    def get_distances(self, ref_object) -> NDArray:
        pass

    @abstractmethod
    def csv(self, variable_names: bool, hash_column: bool,
            attribute_column: bool) -> str:
        """Comma-separated table representation the object array

        Args:
            variable_names: if True first line include the variable names
            hash_column: if True hash will be included as first column
            attribute_column: if True attributes will be included as
                    last column
        Returns:
            CSV representation

        """
        pass

    @abstractmethod
    def mod_round_values(self, decimals: int = 0, int_type: type = np.int16) -> None:
        """Round all values

        Args:
            decimals: number of decimal places
            int_type: numpy int type (default=numpy.int16)
        """
        pass

    def get_distances_matrix(self, between_positions: bool = False) -> NDArray:
        """between position ignores the dot size"""
        if between_positions:
            return spatial.distance.cdist(self._xy, self._xy)
        # matrix with all distance between all points
        dist = np.asarray([self.get_distances(d) for d in self.iter_objects()])
        return dist

    def get_overlaps(self) -> Tuple[np.ndarray, np.ndarray]:
        """return pairs of indices of overlapping of objects and an array of the
        amount of overlap
        takes into account min_dist_between property

        """
        dist = misc.triu_nan(self.get_distances_matrix(between_positions=False),
                             k=1)
        overlap = np.where(dist < self.min_dist_between)
        return np.asarray(overlap).T, np.abs(dist[overlap])

    def get_center_of_mass(self) -> NDArray:
        weighted_sum = np.sum(self._xy * self.perimeter[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self.perimeter)

    def mod_center_array_mass(self) -> None:
        self._xy = self._xy - self.get_center_of_mass()
        self._properties.reset()

    def mod_center_field_area(self) -> None:
        cxy = geometry.center_of_positions(self.properties.convex_hull.xy)
        self._xy = self._xy - cxy
        self._properties.reset()

    def get_free_positionNN(self,
                          ref_object: Union[Dot, Rectangle, Point],
                          in_neighborhood: bool = False,
                          allow_overlapping: bool = False,
                          inside_convex_hull: bool = False,
                          occupied_space=None) -> Union[Dot, Rectangle, Point]:
        """returns the copy of object of at a random free position

        raise exception if not found
        occupied space: see generator generate
        """ #TODO check for point array (point array check is anyway required)

        if isinstance(ref_object, Dot):
            object_size = ref_object.diameter / 2.0
        elif isinstance(ref_object, Rectangle):
            object_size = max(ref_object.size)
        elif isinstance(ref_object, Point):
            object_size = 0
        else:
            raise NotImplementedError("Not implemented for {}".format(
                type(ref_object).__name__))
        if occupied_space is not None and \
                not isinstance(occupied_space, ABCObjectArray):
            raise TypeError("Occupied_space has to be a Dot or Rectangle Array or None.")

        area_rad = self.target_area_radius - self.min_dist_area_boarder - object_size
        rtn_object = deepcopy(ref_object)

        random_walk = BrownianMotion(ref_object.xy, delta=2)

        cnt = 0
        while True:
            if cnt > constants.MAX_ITERATIONS:
                raise NoSolutionError(u"Can't find a free position")
            cnt += 1

            if in_neighborhood:
                rtn_object.xy = random_walk.next()
            else:
                rtn_object.xy = rng.generator.random(size=2) * 2 * area_rad - area_rad

            # is outside area
            if isinstance(ref_object, (Dot, Point)):
                is_outside = area_rad <= rtn_object.polar_radius
            else:
                # Rect: check if one edge is outside
                assert isinstance(rtn_object, Rectangle) # should be alway
                is_outside = False
                for edge in rtn_object.iter_edges():
                    if edge.polar_radius >= area_rad:
                        is_outside = True
                        break

            # check convex hull is not already outside
            if not is_outside and inside_convex_hull:
                # use only those that do not change the convex hull
                tmp_array = self.copy(deep_copy=True)
                tmp_array.add([rtn_object])  # type: ignore
                is_outside = tmp_array.properties.convex_hull != \
                             self.properties.convex_hull
            if is_outside:
                if in_neighborhood:
                    random_walk.step_back()
                continue

            if not allow_overlapping:
                # find overlap (and lower minimum distance)
                dist = self.get_distances(rtn_object)
                if isinstance(occupied_space, ABCObjectArray):
                    dist = np.append(dist,
                                     occupied_space.get_distances(rtn_object))
                if sum(dist < self.min_dist_between) > 0:  # at least one is overlapping
                    continue
            return rtn_object

    def mod_shuffle_positions(self, allow_overlapping: bool = False) -> None:
        """might raise an exception"""
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        all_objects = list(self.iter_objects())
        self.clear()
        for obj in all_objects:
            new = self.get_free_positionNN(obj, in_neighborhood=False,
                            allow_overlapping=allow_overlapping)
            self.add([new])  # type: ignore

    def get_outlier(self) -> Tuple[NDArray, NDArray]:
        """returns indices of object that stand out and array with the size
        of outstanding
        """

        xy = self.properties.convex_hull.xy
        sizes_outlying = np.hypot(xy[:, 0], xy[:, 1]) - \
                         (self.target_area_radius - self.min_dist_area_boarder)
        idx = sizes_outlying > 0

        return (self.properties.convex_hull.object_indices[idx],
                    sizes_outlying[idx])

    def get_number_deviant(self, change_numerosity: int,
                           preserve_field_area: bool = False) -> PointArray:
        """number deviant
        """
        object_array = self.copy()
        new_num = self.properties.numerosity + change_numerosity
        self.properties.fit_numerosity(value=new_num,
                                       keep_convex_hull=preserve_field_area)
        return object_array

    def mod_remove_overlaps(self, keep_convex_hull: bool = False,
                            strict: bool = False) -> bool:
        """
        Returns
        Parameters
        ----------
        keep_convex_hull
        strict

        Returns
        -------
        rtn: boolean
            True, if field area has not been changed (in case strict=False)

        Notes
        -----

        TODO describe different algorithm for keep and not keep CH
        """

        warning_info = "Can't keep field area constant."
        old_fa = self.properties.field_area

        if not keep_convex_hull:
            # touch convex hull objects
            ids = list(range(len(self._xy)))
            for x in self.properties.convex_hull.object_indices_unique:
                self.mod_move_object(x, 0, (0, 0), push_other=True)
                ids.remove(x)
            # touch remaining ids
            for x in ids:
                # touch each object and push other
                self.mod_move_object(x, 0, (0, 0), push_other=True)

        else:
            overlaps = self.get_overlaps()[0]
            ch_idx = self.properties.convex_hull.object_indices_unique

            while len(overlaps):
                # do not replace convex hull objects and try to
                # take that one not on convex hull or raise error/warning
                if overlaps[0, 0] not in ch_idx:
                    idx = overlaps[0, 0]
                elif overlaps[0, 1] not in ch_idx:
                    idx = overlaps[0, 1]
                elif not strict:
                    # warning
                    idx = overlaps[0, 0]
                else:
                    raise NoSolutionError(warning_info)

                obj = next(self.iter_objects(idx))
                self.delete(idx)

                # search new pos: fist inside convex hull later outside
                found = None
                for inside_convex_hull in (True, False):
                    if strict and not inside_convex_hull:
                        continue
                    for in_neighborhood in (True, False):
                        if found is None:
                            try:
                                found = self.get_free_positionNN(
                                    ref_object=obj,
                                    in_neighborhood=in_neighborhood,
                                    inside_convex_hull=inside_convex_hull)
                            except NoSolutionError:
                                found = None

                if found is None:
                    self.add([obj])
                    raise NoSolutionError("Can't find a solution for remove overlap")
                else:
                    self.add([found])
                    overlaps = self.get_overlaps()[0]

            self.properties.reset()

        # check convex hull change
        new_ch = self.properties.field_area
        if keep_convex_hull and old_fa != new_ch:
            if strict:
                raise NoSolutionError(warning_info)
            else:
                print("Warning: " + warning_info)

        return old_fa == new_ch

    def mod_replace(self, xy: Sequence[float]) -> None:
        """Replace the object array

        Args:
          xy: replacement at the x and y axis
        """
        self._xy = self._xy + np.asarray(xy)
        self._properties.reset()

    def mod_move_object(self, object_id: int,
                        distance: float,
                        direction: Union[float, Sequence[float]],
                        push_other: bool = False) -> None:
        """Move a single object by a particular distance. The function uses
        direction and direction

        Args:
            object_id:
            distance:
            direction: polar angle (float) or cartesian coordinate (tuple of two floats)
                       to indicate the direction
            push_other: replace other objects, if required (optional, default=False)
        """
        try:
            ang = float(direction) # type: ignore
        except (TypeError, ValueError):
            try:
                ang = Coordinate(x=direction[0], y=direction[1]) # type: ignore
            except IndexError:
                ang = None
        if ang is None:
            raise TypeError("Direction has to be float or a 2D Coordinate "
                            "or a tuple/list of two elements.")

        obj = next(self.iter_objects(indices=object_id))

        movement = Coordinate()
        if isinstance(ang, float):
            movement.polar = (distance, ang)
        else:
            # "ang" is not an angle it a point
            movement = ang - obj
            movement.polar_radius = distance

        self._xy[object_id, :] = self._xy[object_id, :] + movement.xy

        if push_other:
            # push overlapping object
            obj = next(self.iter_objects(indices=object_id))
            dist = self.get_distances(obj)
            for other_id in np.flatnonzero(dist < 0):
                if other_id != object_id:
                    movement.xy = self._xy[other_id, :] - self._xy[object_id, :]
                    self.mod_move_object(other_id,
                                         direction=movement.polar_angle,
                                         distance=abs(dist[other_id]) + self.min_dist_between,
                                         push_other=True)

        self.properties.reset()

    def mod_squeeze_to_area(self, push_other: bool = True) -> None:
        """Squeeze in target area to remove all standouts

        Args:
            push_other: replace other object, if required to avoid overlaps
        """
        cnt = 0
        while True:
            cnt += 1
            if cnt > constants.MAX_ITERATIONS:
                raise NoSolutionError("Can't find a solution for squeezing")

            idx, size = self.get_outlier()
            if len(idx) == 0:
                return
            for object_id, size in zip(idx, size):
                self.mod_move_object(object_id, distance=size,
                                     direction=(0, 0),
                                     push_other=push_other)

    def get_split_arrays(self) -> List[PointArray]:
        """returns a list of _arrays
        each array contains all dots of with particular colour"""
        att = self._attributes
        att[np.where(att==None)] = "None"

        rtn = []
        for c in np.unique(att):
            if c is not None:
                da = self.copy(indices=0, deep_copy=False)  # fast. shallow copy with just one object
                da.clear()
                da.add(self.find_objects(attribute=c))
                rtn.append(da)
        return rtn

    def mod_realign(self, keep_convex_hull=False, strict=True) -> Tuple[bool, bool]:
        """

        Parameters
        ----------
        keep_convex_hull
        strict

        Returns
        -------
        field_area_unchanged, no_outlier

        """
        convex_hull_unchanged = self.mod_remove_overlaps(keep_convex_hull=keep_convex_hull,
                                                         strict=strict)
        self.mod_center_field_area()

        has_no_outlier = len(self.get_outlier()[0]) == 0
        if not has_no_outlier:
            if keep_convex_hull:
                warning_info = "Can't keep field area constant."
                if strict:
                    raise NoSolutionError(warning_info)
                else:
                    print("Warning: " + warning_info)

            self.mod_squeeze_to_area(push_other=True)
            has_no_outlier = True
            convex_hull_unchanged = False

        return convex_hull_unchanged, has_no_outlier

# TODO  everywhere: file header doc and author information
