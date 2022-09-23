"""

"""
from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

import json
from abc import ABCMeta, abstractmethod
from copy import deepcopy
from hashlib import md5
from typing import Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import spatial

from .._lib import geometry, np_tools, rng
from .._lib.coordinate import Coordinate
from .._lib.exception import NoSolutionError
from .._lib.misc import dict_to_text
from .._lib.typing import IntOVector
from .._shapes import ShapeType
from .._shapes.dot import Dot
from .._shapes.rectangle import Rectangle
from .._shapes.dot_list import DotList
from .._shapes.rectangle_list import RectangleList
from .properties import ArrayProperties
from .target_area import TargetArea
from .tools import BrownianMotion
from .. import constants


class ABCObjectArray(TargetArea, metaclass=ABCMeta):
    """Abstract array class

    All object arrays are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(self,
                 target_area: ShapeType,
                 min_distance_between_objects: Optional[float],
                 min_distance_area_boarder: Optional[float],
                 object_shapes: Union[DotList, RectangleList],
                 attributes: Optional[ArrayLike] = None) -> None:
        # also as parent  class for implementation of dot and rect arrays

        super().__init__(target_area=target_area,
                         min_distance_between_objects=min_distance_between_objects,
                         min_distance_area_boarder=min_distance_area_boarder)

        self._objs = object_shapes
        self._attributes = np.array([])
        self._set_missing_attributes(attributes)
        self._properties = ArrayProperties(self)

    @ property
    @ abstractmethod
    def surface_areas(self):
        """Vector with the surface areas of all objects"""

    @ property
    @ abstractmethod
    def perimeter(self) -> NDArray:
        """Vector with the perimeter of all objects"""

    @ abstractmethod
    def to_dict(self) -> dict:
        """Convert array to dictionary
        """
        d = super().to_dict()
        d.update({"xy": self._objs.xy.tolist()})
        if len(self._attributes) > 0:
            att = self._attributes.flatten()
            if np.all(att == att[0]):  # all equal
                d.update({"attributes": self._attributes[0]})
            else:
                d.update({"attributes": self._attributes.tolist()})
        return d

    @ abstractmethod
    def dataframe_dict(self):
        """Returns dict that can be used to create Pandas dataframe or Arrow Table

        Examples
        --------
        >>> df_dict = stimulus.dataframe_dict()
        >>> df = pandas.DataFame(df_dict) # Pandas dataframe

        >>> table = pyarrow.Table.from_pydict(df_dict) # Arrow Table
        """

    @ staticmethod
    @ abstractmethod
    def from_dict(the_dict):
        """ """

    @ abstractmethod
    def copy(self, indices: Union[IntOVector, None] = None,
             deep_copy=True) -> ABCObjectArray:
        """ """

    @ abstractmethod
    def iter_objects(self, indices: Optional[IntOVector] = None) \
            -> Iterator[ShapeType]:
        """iterate over all or a part of the objects

        Parameters
        ----------
        indices: int or iterable of integer

        Notes
        -----
        To iterate all object you might all use the class iterator __iter__:
        >>> for obj in my_array:
        >>>    print(obj)
        """

    @ abstractmethod
    def add(self, obj):
        """ """

    @ abstractmethod
    def find_objects(self, size=None, attribute=None):
        """ """

    @ abstractmethod
    def get_distances(self, ref_object) -> NDArray:
        """get_distances """

    @ abstractmethod
    def csv(self, variable_names: bool, hash_column: bool,
            attribute_column: bool):
        """Comma-separated table representation the object array

        Args:
            variable_names: if True first line include the variable names
            hash_column: if True hash will be included as first column
            attribute_column: if True attributes will be included as
                    last column
        Returns:
            CSV representation

        """

    @ abstractmethod
    def mod_round_values(self, decimals: int = 0,
                         int_type: type = np.int16) -> None:
        """Round all values

        Args:
            decimals: number of decimal places
            int_type: numpy int type (default=numpy.int16)
        """

    def __str__(self) -> str:
        d = TargetArea.to_dict(self)  # super: omit objects
        prop = self.properties.as_text(extended_format=True)
        return dict_to_text(d, col_a=30, col_b=7) + "\n " + prop[1:]

    @ property
    def xy(self) -> NDArray:
        """Numpy array of the object locations

        The two dimensional array (shape=[2, `n`]) represents the locations of the center of
        the `n` objects in this array
        """
        return self._objs.xy

    @ property
    def attributes(self) -> NDArray:
        """Numpy vector of the object attributes
        """
        return self._attributes

    def set_attributes(self, attributes: ArrayLike) -> None:
        """Set all attributes

        Args:
            attributes: attribute (string) or list of attributes
        """
        if isinstance(attributes, (np.ndarray, list, tuple)):
            if len(attributes) != self._properties.numerosity:
                raise ValueError("Length of attribute list does not adapt the " +
                                 "size of the dot array.")
            self._attributes = np.array(attributes)
        else:
            self._attributes = np.array(
                [attributes] * self._properties.numerosity)

    def _set_missing_attributes(self, attributes: Optional[ArrayLike]):
        """appends attributes. If no list or np.array set all missing attributes
        """
        n_missing_attributes = self._objs.xy.shape[0] - \
            self._attributes.shape[0]
        if not isinstance(attributes, (np.ndarray, tuple, list)):
            attributes = np.array([attributes] * n_missing_attributes)
        elif len(attributes) != n_missing_attributes:
            raise RuntimeError("Badly shaped data: attributes have not " +
                               "the same length as the xy coordinates")

        self._attributes = np.append(self._attributes, attributes)

#    def _append_xy_attribute(self, xy: ArrayLike,
#                             attributes: Optional[ArrayLike] = None) -> None:
#        """appends xy and attributes"""
#        self._objs.append(xy)
#        self._set_missing_attributes(attributes)
#        self._properties.reset()

    def clear(self) -> None:
        self._objs.clear()
        self._attributes = np.array([])
        self._properties.reset()

    def delete(self, index: IntOVector) -> None:
        """Delete object

        Args:
            index:

        Returns:

        """
        self._objs.delete(index)
        self._attributes = np.delete(self._attributes, index)
        self._properties.reset()

    @ property
    def properties(self) -> ArrayProperties:
        """Properties of the object array.

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

    @ property
    def hash(self) -> str:
        """Hash (MD5 hash) of the array

        The hash can be used as an unique identifier of the object array.

        Notes:
            Hashing is based on the byte representations of the positions, perimeter
            and attributes.
        """
        m = md5()
        m.update(
            self._objs.xy.tobytes())  # to_byte required: https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        try:
            m.update(self.perimeter.tobytes())
        except AttributeError:
            pass
        m.update(self.attributes.tobytes())
        return m.hexdigest()

    def to_json(self, file_path: Optional[str] = None,
                indent: Optional[int] = None,
                include_hash: bool = False) -> Union[str, None]:
        """ """

        d = self.to_dict()
        if include_hash:
            d.update({"hash": self.hash})
        j = json.dumps(d, indent=indent)
        if file_path is None:
            return j
        else:
            with open(file_path, 'w', encoding="utf-8") as fl:
                fl.write(j)

    def get_objects(self, indices: Optional[IntOVector] = None) \
            -> Sequence[ShapeType]:
        """ """
        return list(self.iter_objects(indices=indices))

    def get_object(self, index: int) -> Union[None, Dot, Rectangle]:
        """ """
        if isinstance(index, int):
            return next(self.iter_objects(indices=index))
        else:
            raise ValueError("Index must be a integer not a {}. ".format(
                type(index).__name__) + "To handle multiple indices use 'get_objects'. ")

    def join(self, object_array) -> None:
        """add another object _arrays"""
        self.add(object_array.iter_objects())

    def get_distances_matrix(self, between_positions: bool = False) -> NDArray:
        """between position ignores the dot size"""
        if between_positions:
            return spatial.distance.cdist(self._objs.xy, self._objs.xy)
        # matrix with all distance between all objects
        dist = np.asarray([self.get_distances(d) for d in self.iter_objects()])
        return dist

    def get_overlaps(self) -> Tuple[NDArray, NDArray]:
        """return pairs of indices of overlapping of objects and an array of the
        amount of overlap
        takes into account min_distance_between_objects property
        """
        dist = np_tools.triu_nan(self.get_distances_matrix(between_positions=False),
                                 k=1)
        overlap = np.where(dist < self.min_distance_between_objects)
        return np.asarray(overlap).T, np.abs(dist[overlap])

    def get_center_of_mass(self) -> NDArray:
        weighted_sum = np.sum(
            self._objs.xy * self.perimeter[:, np.newaxis], axis=0)
        return weighted_sum / np.sum(self.perimeter)

    def get_center_of_field_area(self) -> NDArray:
        """Center of all objects

        Returns:
            Coordinate of center position
        """
        return geometry.center_of_coordinates(self.properties.convex_hull.xy)

    def mod_center_array_mass(self) -> None:
        self._objs.xy = self._objs.xy - self.get_center_of_mass()
        self._properties.reset()

    def mod_center_field_area(self) -> None:
        cxy = geometry.center_of_coordinates(
            self.properties.convex_hull.xy)
        self._objs.xy = self._objs.xy - cxy  # type: ignore
        self._properties.reset()

    def _search_area_parameter(self):
        """helper function that returns
            the area of possible object location and

        It takes into account to minimum distance to area boarder
        """
        if isinstance(self.target_area, Dot):
            tmp = self.target_area.diameter/2.0 \
                - self.min_distance_area_boarder
            search_area = Dot(xy=(0, 0), diameter=tmp*2)
            half_search_area_size = np.array([tmp, tmp])

        elif isinstance(self.target_area, Rectangle):
            tmp = (self.target_area.width - 2 * self.min_distance_area_boarder,
                   self.target_area.height - 2 * self.min_distance_area_boarder)
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
            raise NotImplementedError("Not implemented for {}".format(
                type(ref_object).__name__))

        if occupied_space is not None and \
                not isinstance(occupied_space, ABCObjectArray):
            raise TypeError(
                "Occupied_space has to be a Dot or Rectangle Array or None.")

        search_area, half_search_area_size = self._search_area_parameter()
        rtn_object = deepcopy(ref_object)
        random_walk = BrownianMotion(ref_object.xy, delta=2)

        cnt = 0
        while True:
            if cnt > constants.MAX_ITERATIONS:
                raise NoSolutionError("Can't find a free position: "
                                      + f"Current n={len(self.xy)}")
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
                tmp_array.add([rtn_object])  # type: ignore
                is_inside = tmp_array.properties.convex_hull == \
                    self.properties.convex_hull

            if not is_inside:
                if in_neighborhood:
                    random_walk.step_back()
                continue  # try another position

            if not allow_overlapping:
                # find overlap (and lower minimum distance)
                dist = self.get_distances(rtn_object)
                if isinstance(occupied_space, ABCObjectArray):
                    dist = np.append(dist,
                                     occupied_space.get_distances(rtn_object))
                if sum(dist < self.min_distance_between_objects) > 0:  # at least one is overlapping
                    continue  # try another position
            return rtn_object

    def mod_shuffle_positions(self, allow_overlapping: bool = False) -> None:
        """might raise an exception"""
        # find new position for each dot
        # mixes always all position (ignores dot limitation)

        all_objects = list(self.iter_objects())
        self.clear()
        for obj in all_objects:
            new = self.get_free_position(obj, in_neighborhood=False,
                                         allow_overlapping=allow_overlapping)
            self.add([new])  # type: ignore

    def get_outlier(self) -> List[int]:
        """Indices of object that are not fully inside the target area"""
        # TODO test me
        indices = []
        search_area, _ = self._search_area_parameter()

        for cnt, obj in enumerate(self.iter_objects()):
            if not obj.is_inside(search_area):
                indices.append(cnt)
        return indices

    def get_number_deviant(self, change_numerosity: int,
                           preserve_field_area: bool = False) -> ABCObjectArray:
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
            ids = list(range(len(self._objs.xy)))
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
                                found = self.get_free_position(
                                    ref_object=obj,
                                    in_neighborhood=in_neighborhood,
                                    inside_convex_hull=inside_convex_hull)
                            except NoSolutionError:
                                found = None

                if found is None:
                    self.add([obj])  # put back
                    raise NoSolutionError(
                        "Can't find a solution for remove overlap")
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

    def mod_replace(self, xy: ArrayLike) -> None:
        """Replace the object array

        Args:
          xy: replacement at the x and y axis
        """
        self._objs.xy = self._objs.xy + np.atleast_2d(xy)
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

        for id_, obj in zip(object_ids, self.iter_objects(indices=object_ids)):

            if isinstance(ang, float):
                movement.polar = (distance, ang)
            else:
                # "ang" is not an angle it an object
                movement.xy = ang.xy - obj.xy
                movement.rho = distance

            self._objs.xy[id_,
                          :] = self._objs.xy[id_, :] + movement.xy

            if push_other:
                # push overlapping object
                obj.xy += movement.xy
                dist = self.get_distances(obj)
                for other_id in np.flatnonzero(dist < self.min_distance_between_objects):
                    if other_id != id_:
                        movement.xy = self._objs.xy[other_id, :] \
                            - self._objs.xy[id_, :]
                        self.mod_move_object(other_id,
                                             direction=movement.theta,
                                             distance=abs(dist[other_id])
                                             + self.min_distance_between_objects,
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
            for object_id in idx:
                obj = self.get_object(index=object_id)
                if isinstance(self.target_area, Dot):
                    mv_dist = obj.distance(self.target_area)  # type: ignore
                    print(mv_dist)
                    self.mod_move_object(object_id,
                                         distance=abs(mv_dist),
                                         direction=(0, 0),
                                         push_other=push_other)
                else:
                    # target area is rect
                    # should be always this case
                    assert isinstance(self.target_area, Rectangle)
                    raise NotImplementedError()  # FIXME

    def get_split_arrays(self) -> List[ABCObjectArray]:
        """returns a list of _arrays
        each array contains all dots of with particular colour"""
        att = self._attributes
        att[np.where(att == None)] = "None"

        rtn = []
        for c in np.unique(att):
            if c is not None:
                # fast. shallow copy with just one object
                da = self.copy(indices=0, deep_copy=False)
                da.clear()
                da.add(self.find_objects(attribute=c))  # type: ignore
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

# TODO  everywhere: file header doc and author information
