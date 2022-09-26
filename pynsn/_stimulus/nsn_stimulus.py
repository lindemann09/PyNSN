"""

"""
from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import json
from copy import deepcopy
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import spatial

from .. import _stimulus, constants
from .._lib import np_tools, rng
from .._lib.exceptions import NoSolutionError
from .._lib.misc import dict_to_text
from .._lib.typing import IntOVector
from .._object_arrays import (BaseDotArray, BaseRectangleArray, DotArray,
                              RectangleArray)
from .._shapes import Coordinate, Dot, Rectangle, ShapeType
from ..spatial_relations import _spatial_relations as sprel
from .target_area import TargetArea
from .tools import BrownianMotion

DOT_ARRAY = "Dot Array"
RECTANGLE_ARRAY = "Rectangle Array"


class NSNStimulus(TargetArea):
    """Non-Symbolic Number Stimulus

    NSN-Stimulus are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(self,
                 object_array: Union[DotArray, RectangleArray],
                 target_area_shape: ShapeType,
                 min_dist_between_objects: Optional[float] = None,
                 min_dist_area_edge: Optional[float] = None) -> None:
        # also as parent  class for implementation of dot and nsn stimuli

        super().__init__(target_area_shape=target_area_shape,
                         min_dist_between_objects=min_dist_between_objects,
                         min_dist_area_edge=min_dist_area_edge)
        if not isinstance(object_array, (DotArray, RectangleArray)):
            raise ValueError("object_array has to be a DotArray "
                             "or RectangleArray, but not "
                             f"{type(object_array).__name__}")
        self._objects = object_array
        self._properties = _stimulus.ArrayProperties(self)

    @property
    def objects(self) -> Union[DotArray, RectangleArray]:
        """The objects in the array"""
        return self._objects

    def type(self) -> str:
        if isinstance(self.objects, DotArray):
            return DOT_ARRAY
        elif isinstance(self.objects, RectangleArray):
            return RECTANGLE_ARRAY
        else:
            raise ValueError("object_array has to be a DotArray "
                             "or RectangleArray, but not "
                             f"{type(self.objects).__name__}")

    def to_dict(self) -> dict:
        """Convert stimulus to dictionary
        """
        rtn = {"type": self.type()}
        rtn.update(super().to_dict())
        rtn.update(self.objects.to_dict())
        return rtn

    @classmethod
    def from_dict(cls, the_dict: dict) -> NSNStimulus:
        """read nsn stimulus from dict"""

        oa_type = the_dict["type"]
        if oa_type == DOT_ARRAY:
            object_array = DotArray.from_dict(the_dict)
        elif oa_type == RECTANGLE_ARRAY:
            object_array = RectangleArray.from_dict(the_dict)
        else:
            raise RuntimeError(f"Unknown nsn stimulus type {oa_type}")
        target_area = TargetArea.from_dict(the_dict)

        return cls(
            object_array=object_array,
            target_area_shape=target_area.target_area_shape,
            min_dist_between_objects=target_area.min_dist_between_objects,
            min_dist_area_edge=target_area.min_dist_area_edge)

    @classmethod
    def from_json(cls, filename_or_json_str: str) -> NSNStimulus:
        """Loading stimulus from json file or json string

        Args:
            filename: the file name

        Returns:
            nsn stimulus
        """
        fl = Path(filename_or_json_str)
        if fl.is_file():
            json_str = fl.read_text(encoding="utf-8")
        else:
            json_str = filename_or_json_str
        return cls.from_dict(json.loads(json_str))

    def copy(self, indices: Union[int, Sequence[int], None] = None,
             deep_copy: bool = True) -> NSNStimulus:
        """returns a (deep) copy of the nsn stimulus.

        It allows to copy a subset of objects.

        """
        if deep_copy:
            target_area_shape = deepcopy(self.target_area_shape)
        else:
            target_area_shape = self.target_area_shape
        return NSNStimulus(
            object_array=self._objects.copy(
                indices=indices, deep_copy=deep_copy),
            target_area_shape=target_area_shape,
            min_dist_between_objects=self.min_dist_between_objects,
            min_dist_area_edge=self.min_dist_area_edge)

    def __str__(self) -> str:
        d = TargetArea.to_dict(self)  # super: omit objects
        prop = self.properties.as_text(extended_format=True)

        rtn = f"-{self.type():-^38}\n"
        rtn += f" {self._objects.hash():^38}\n "
        return rtn + dict_to_text(d)[1:] + f"\n {'Properties':^38}\n " + prop[1:]

    def clear(self) -> None:
        self._objects.clear()
        self._properties.reset()

    def delete(self, index: IntOVector) -> None:
        """Delete object

        Args:
            index:

        Returns:

        """
        self._objects.delete(index)
        self._properties.reset()

    @property
    def properties(self) -> _stimulus.ArrayProperties:
        """Properties of the nsn stimulus.

        ``ArrayProperties`` represents and handles (fitting, scaling) visual
        properties of the object like

        * numerosity
        * average_dot_diameter/average_rectangle_size
        * total_surface_area
        * average_surface_area
        * total_perimeter
        * average_perimeter
        * field_area
        * field_area_positions
        * sparsity
        * log_spacing
        * log_size
        * coverage
        """
        return self._properties

    def to_json(self, file_path: Optional[str] = None,
                indent: Optional[int] = None,
                include_hash: bool = False) -> Union[str, None]:
        """ """

        d = self.to_dict()
        if include_hash:
            d.update({"hash": self._objects.hash()})

        json_str = json.dumps(d, indent=indent)
        if file_path is not None:
            with open(file_path, 'w', encoding="utf-8") as fl:
                fl.write(json_str)
        return json_str

    def join(self, nsn_stimulus) -> None:
        """joins with another dot array """
        self._objects.join(nsn_stimulus.objects)
        self._properties.reset()

    def get_distances_matrix(self, between_positions: bool = False) -> NDArray:
        """between position ignores the dot size"""
        if between_positions:
            return spatial.distance.cdist(self._objects.xy, self._objects.xy)
        # matrix with all distance between all objects
        dist = np.asarray([self.get_distances(d)
                           for d in self._objects.iter()])
        return dist

    def get_overlaps(self) -> Tuple[NDArray, NDArray]:
        """return pairs of indices of overlapping of objects and an array of the
        amount of overlap
        takes into account min_dist_between_objects property
        """
        dist = np_tools.triu_nan(self.get_distances_matrix(between_positions=False),
                                 k=1)
        overlap = np.where(dist < self.min_dist_between_objects)
        return np.asarray(overlap).T, np.abs(dist[overlap])

    def mod_center_array_mass(self) -> None:
        self._objects.xy = self._objects.xy - self._objects.center_of_mass
        self._properties.reset()

    def mod_center_field_area(self) -> None:
        cxy = self.properties.convex_hull.center
        self._objects.xy = self._objects.xy - cxy  # type: ignore
        self._properties.reset()

    def _search_area_parameter(self):
        """helper function that returns
            the area of possible object location and

        It takes into account to minimum distance to area boarder
        """
        if isinstance(self.target_area_shape, Dot):
            tmp = self.target_area_shape.diameter/2.0 \
                - self.min_dist_area_edge
            search_area = Dot(xy=(0, 0), diameter=tmp*2)
            half_search_area_size = np.array([tmp, tmp])

        elif isinstance(self.target_area_shape, Rectangle):
            tmp = (self.target_area_shape.width - 2 * self.min_dist_area_edge,
                   self.target_area_shape.height - 2 * self.min_dist_area_edge)
            search_area = Rectangle(xy=(0, 0), size=tmp)
            half_search_area_size = np.asarray(search_area.size) / 2.0
        else:
            raise NotImplementedError()

        return search_area, half_search_area_size

    def get_free_position(self,
                          ref_object: ShapeType,
                          in_neighborhood: bool = False,
                          allow_overlapping: bool = False,
                          inside_convex_hull: bool = False,
                          occupied_space=None) -> ShapeType:
        """returns the copy of object of at a random free position

        raise exception if not found
        occupied space: see generator generate
        """

        if not isinstance(ref_object, (Dot, Rectangle)):
            raise NotImplementedError("Not implemented for "
                                      f"{type(ref_object).__name__}")

        if occupied_space is not None and \
                not isinstance(occupied_space, NSNStimulus):
            raise TypeError(
                "Occupied_space has to be a Dot or nsn stimulus or None.")

        search_area, half_search_area_size = self._search_area_parameter()
        rtn_object = deepcopy(ref_object)
        random_walk = BrownianMotion(ref_object.xy, delta=2)

        cnt = 0
        while True:
            if cnt > constants.MAX_ITERATIONS:
                raise NoSolutionError("Can't find a free position: "
                                      + f"Current n={len(self._objects.xy)}")
            cnt += 1

            if in_neighborhood:
                rtn_object.xy = random_walk.next()
            else:
                rtn_object.xy = rng.generator.random(size=2) * 2 \
                    * half_search_area_size - half_search_area_size

            # is outside area
            is_inside = rtn_object.is_inside(search_area)
            # check convex hull is not already outside
            if is_inside and inside_convex_hull:
                # use only those that do not change the convex hull
                tmp_array = self.copy(deep_copy=True)
                tmp_array.objects.add([rtn_object])  # type: ignore
                is_inside = tmp_array.properties.convex_hull == \
                    self.properties.convex_hull

            if not is_inside:
                if in_neighborhood:
                    random_walk.step_back()
                continue  # try another position

            if not allow_overlapping:
                # find overlap (and lower minimum distance)
                dist = self.get_distances(rtn_object)
                if isinstance(occupied_space, NSNStimulus):
                    dist = np.append(dist,
                                     occupied_space.get_distances(rtn_object))
                if sum(dist < self.min_dist_between_objects) > 0:  # at least one is overlapping
                    continue  # try another position
            return rtn_object

    def mod_shuffle_positions(self, allow_overlapping: bool = False) -> None:
        """might raise an exception"""
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        all_objects = list(self._objects.iter())
        self.clear()
        for obj in all_objects:
            new = self.get_free_position(obj, in_neighborhood=False,
                                         allow_overlapping=allow_overlapping)
            self.objects.add([new])  # type: ignore

    def get_outlier(self) -> List[int]:
        """Indices of object that are not fully inside the target area"""
        # FIXME make numpy search test me
        indices = []
        search_area, _ = self._search_area_parameter()

        for cnt, obj in enumerate(self._objects.iter()):
            if not obj.is_inside(search_area):
                indices.append(cnt)
        return indices

    def get_number_deviant(self, change_numerosity: int,
                           preserve_field_area: bool = False) -> NSNStimulus:
        """number deviant
        """
        nsn_stimulus = self.copy()
        new_num = self.properties.numerosity + change_numerosity
        self.properties.fit_numerosity(value=new_num,
                                       keep_convex_hull=preserve_field_area)
        return nsn_stimulus

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
        # FIXME not tested
        warning_info = "Can't keep field area constant."
        old_fa = self.properties.field_area

        if not keep_convex_hull:
            # touch convex hull objects
            ids = list(range(len(self._objects.xy)))
            for x in self.properties.convex_hull.object_indices:
                self.mod_move_object(x, 0, (0, 0), push_other=True)
                ids.remove(x)
            # touch remaining ids
            for x in ids:
                # touch each object and push other
                self.mod_move_object(x, 0, (0, 0), push_other=True)

        else:
            overlaps = self.get_overlaps()[0]
            ch_idx = self.properties.convex_hull.object_indices

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

                obj = next(self._objects.iter(idx))
                self.delete(idx)

                # search new pos: fist inside convex hull later outside
                found = None
                for inside_convex_hull in (True, False):
                    if strict and not inside_convex_hull:
                        continue
                    for in_neighborhood in (True, False):
                        if found is None:
                            try:
                                found = self.get_free_position(
                                    ref_object=obj,
                                    in_neighborhood=in_neighborhood,
                                    inside_convex_hull=inside_convex_hull)
                            except NoSolutionError:
                                found = None

                if found is None:
                    self.objects.add([obj])  # type: ignore # put back
                    raise NoSolutionError(
                        "Can't find a solution for remove overlap")
                else:
                    self.objects.add([found])
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

    def mod_replace(self, xy: ArrayLike) -> None:
        """Replace the nsn stimulus

        Args:
          xy: replacement at the x and y axis
        """
        self._objects.xy = self._objects.xy + np.atleast_2d(xy)
        self._properties.reset()

    def mod_move_object(self, object_ids: IntOVector,
                        distance: float,
                        direction: Union[float, ArrayLike],
                        push_other: bool = False) -> None:
        """Move a single object by a particular distance. The function uses
        direction and direction

        Args:
            object_id: the object to be moved
            distance: replacement distance
            direction: polar angle (float) or cartesian coordinate (tuple of two floats)
                       toward all objects should be moved. In the last case, the movement
                       direction differs for each object.
            push_other: replace other objects, if required (optional, default=False)
        """
        if isinstance(direction, float):
            ang = direction
        else:
            try:
                ang = Coordinate(xy=direction)
            except IndexError:
                ang = None
        if ang is None:
            raise ValueError("Direction has to be float or a 2D Coordinate "
                             ", thus, an ArrayLike with two elements.")
        if isinstance(object_ids, (int, np.integer)):
            object_ids = [object_ids]
        movement = Coordinate(xy=(0, 0))

        for id_, obj in zip(object_ids, self._objects.iter(indices=object_ids)):

            if isinstance(ang, float):
                movement.polar = (distance, ang)
            else:
                # "ang" is not an angle it an object
                movement.xy = ang.xy - obj.xy
                movement.rho = distance

            self._objects.xy[id_,
                             :] = self._objects.xy[id_, :] + movement.xy

            if push_other:
                # push overlapping object
                obj.xy += movement.xy
                dist = self.get_distances(obj)
                for other_id in np.flatnonzero(dist < self.min_dist_between_objects):
                    if other_id != id_:
                        movement.xy = self._objects.xy[other_id, :] \
                            - self._objects.xy[id_, :]
                        self.mod_move_object(other_id,
                                             direction=movement.theta,
                                             distance=abs(dist[other_id])
                                             + self.min_dist_between_objects,
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

            idx = self.get_outlier()
            if len(idx) == 0:
                return
            for object_id, obj in zip(idx, self.objects.iter(idx)):
                if isinstance(self.target_area_shape, Dot):
                    mv_dist = obj.distance(
                        self.target_area_shape)  # type: ignore
                    print(mv_dist)
                    self.mod_move_object(object_id,
                                         distance=abs(mv_dist),
                                         direction=(0, 0),
                                         push_other=push_other)
                else:
                    # target area is rect
                    # should be always this case
                    assert isinstance(self.target_area_shape, Rectangle)
                    raise NotImplementedError()  # FIXME

    def get_split_arrays(self) -> List[NSNStimulus]:
        """returns a list of nsn stimuli
        each array contains all dots of with particular colour"""
        # FIXME not checked
        att = self._objects.attributes
        att[np.where(att == None)] = "None"

        rtn = []
        for c in np.unique(att):
            if c is not None:
                # fast. shallow copy with just one object
                da = self.copy(indices=0, deep_copy=False)
                da.clear()
                idx = self.objects.find(attribute=c)
                objects = list(self._objects.iter(indices=idx))
                da.objects.add(objects)  # type: ignore
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

        has_no_outlier = len(self.get_outlier()) == 0
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

    def get_distances(self, shape: ShapeType) -> NDArray:
        """Euclidean distances toward a single object

        Negative numbers indicate an overlap of the objects

        Args:
            dot: object towards the distances will be calculated

        Returns: numpy array of distances
        """

        if len(self._objects.xy) == 0:
            return np.array([])

        if isinstance(shape, Dot) and isinstance(self._objects, DotArray):
            rel = sprel.DotDot(dots_a=self._objects,
                               dots_b=BaseDotArray(xy=shape.xy,
                                                   diameter=shape.diameter))
            return rel.distances()

        elif isinstance(shape, Rectangle) and \
                isinstance(self._objects, RectangleArray):
            rel = sprel.RectangleRectangle(
                rectangles_a=self._objects,
                rectangles_b=BaseRectangleArray(xy=shape.xy, sizes=shape.size))
            return rel.xy_distances()

        else:
            raise NotImplementedError()
