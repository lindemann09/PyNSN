"""
Dot Array
"""
from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import deepcopy
from typing import Any, Dict, Iterator, List, Optional, Sequence, Union

import numpy as np
from numpy.typing import NDArray, ArrayLike
from .._lib import np_tools
from .._lib.spatial_relations import DotDot
from .._shapes import Dot,  ShapeType
from .._shapes.np_shapes import NPDots
from .abc_object_array import ABCObjectArray, IntOVector, NPDots
from .tools import make_csv

# pylint: disable=W0237:arguments-renamed


# TODO: How to deal with rounding? Is saving to precises? Suggestion:
#  introduction precision parameter that is used by to_dict and get_csv and
#  hash


class DotArray(ABCObjectArray):
    """Dot array class

    All object arrays are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(self,
                 target_area: ShapeType,
                 min_distance_between_objects: Optional[float] = None,
                 min_distance_area_boarder: Optional[float] = None,
                 xy: Optional[ArrayLike] = None,
                 diameter: Optional[ArrayLike] = None,
                 attributes: Optional[ArrayLike] = None) -> None:

        self._np_shapes = NPDots(xy, diameter)  # for typing
        super().__init__(np_shapes=self._np_shapes,
                         attributes=attributes,
                         target_area=target_area,
                         min_distance_between_objects=min_distance_between_objects,
                         min_distance_area_boarder=min_distance_area_boarder)

    def add(self, obj: Union[Dot, Sequence[Dot]]) -> None:
        """append one dot or list of dots"""
        if isinstance(obj, Dot):
            obj = [obj]
        else:
            obj = list(obj)

        xy = []
        diameter = []
        attributes = []
        for d in obj:
            xy.append(d.xy)
            attributes.append(d.attribute)
            diameter.append(d.diameter)

        self._np_shapes.append(xy, diameter)
        self._set_missing_attributes(attributes)
        self.properties.reset()

    @property
    def diameter(self) -> NDArray:
        return self._np_shapes.diameter

    @property
    def surface_areas(self) -> NDArray:
        """TODO

        """
        # a = pi r**2 = pi d**2 / 4
        return np.pi * (self._np_shapes.diameter ** 2) / 4.0

    @property
    def perimeter(self) -> NDArray:
        """Perimeter for each dot

        """
        return np.pi * self._np_shapes.diameter

    def mod_round_values(self, decimals: int = 0,
                         int_type: type = np.int32) -> None:
        # inherited doc
        if decimals is None:
            return
        self._np_shapes.xy = np_tools.round2(self._np_shapes.xy, decimals=decimals,
                                             int_type=int_type)
        self._np_shapes.diameter = np_tools.round2(self._np_shapes.diameter, decimals=decimals,
                                                   int_type=int_type)

    def to_dict(self) -> dict:
        d = super().to_dict()
        d.update({"diameter": self._np_shapes.diameter.tolist()})
        return d

    def dataframe_dict(self, hash_column: bool = False,
                       attribute_column: bool = True) -> dict:
        # inherited doc
        if hash_column:
            d = {"hash": [self.hash] * len(self._np_shapes.xy)}
        else:
            d = {}
        d.update({"x": self._np_shapes.xy[:, 0].tolist(),
                  "y": self._np_shapes.xy[:, 1].tolist(),
                  "diameter": self._np_shapes.diameter.tolist()})
        if attribute_column:
            d.update({"attributes": self._attributes.tolist()})
        return d

    @staticmethod
    def from_dict(the_dict: Dict[str, Any]) -> DotArray:
        """read Dot collection from dict"""
        rtn = DotArray(target_area=DotArray._target_area_from_dict(the_dict),
                       min_distance_between_objects=the_dict["min_distance_between_objects"],
                       min_distance_area_boarder=the_dict["min_distance_area_boarder"])

        rtn._np_shapes.append(xy=the_dict["xy"], diameter=the_dict["diameter"])
        rtn._set_missing_attributes(attributes=the_dict["attributes"])
        return rtn

    def copy(self, indices: Union[int, Sequence[int], None] = None,
             deep_copy: bool = True) -> DotArray:
        """A (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        Args:
            indices: Arrays with the indices to be be copied. If None, all
                object will be copies

        Returns:
            a copy of the array
        """

        if len(self._np_shapes.xy) == 0:
            return DotArray(target_area=deepcopy(self.target_area),
                            min_distance_between_objects=self.min_distance_between_objects,
                            min_distance_area_boarder=self.min_distance_area_boarder)

        if indices is None:
            indices = list(range(len(self._np_shapes.xy)))

        if deep_copy:
            return DotArray(target_area=deepcopy(self.target_area),
                            min_distance_between_objects=self.min_distance_between_objects,
                            min_distance_area_boarder=self.min_distance_area_boarder,
                            xy=self._np_shapes.xy[indices, :].copy(),
                            diameter=self._np_shapes.diameter[indices].copy(),
                            attributes=self._attributes[indices].copy())
        else:
            return DotArray(target_area=deepcopy(self.target_area),
                            min_distance_between_objects=self.min_distance_between_objects,
                            min_distance_area_boarder=self.min_distance_area_boarder,
                            xy=self._np_shapes.xy[indices, :],
                            diameter=self._np_shapes.diameter[indices],
                            attributes=self._attributes[indices])

    def get_distances(self, dot: Dot) -> NDArray:
        """Euclidean distances toward a single dot

        Negative numbers indicate an overlap of the objects

        Args:
            dot: object towards the distances will be calculated

        Returns: numpy array of distances
        """
        assert isinstance(dot, Dot)
        if len(self._np_shapes.xy) == 0:
            return np.array([])
        else:
            rel = DotDot(a_xy=self._np_shapes.xy,
                         a_diameter=self._np_shapes.diameter,
                         b_xy=dot.xy,
                         b_diameter=np.array([dot.diameter]))
            return rel.distances()

    def iter_objects(self, indices: Optional[IntOVector] = None) -> Iterator[Dot]:
        """iterate over all or a part of the objects

        Parameters
        ----------
        indices: int or interable of integer

        Notes
        -----
        To iterate all object you might all use the class iterator __iter__:
        >>> for obj in my_array:
        >>>    print(obj)
        """

        if isinstance(indices, (int, np.integer)):
            yield Dot(xy=self._np_shapes.xy[indices, :],
                      diameter=self._np_shapes.diameter[indices],
                      attribute=self._attributes[indices])
        else:
            if indices is None:
                data = zip(self._np_shapes.xy,
                           self._np_shapes.diameter,
                           self._attributes)
            else:
                data = zip(self._np_shapes.xy[indices, :],  # type: ignore
                           self._np_shapes.diameter[indices],
                           self._attributes[indices])
            for xy, dia, att in data:
                yield Dot(xy=xy, diameter=dia, attribute=att)

    def find_objects(self, diameter: Optional[float] = None,
                     attribute: Any = None) -> List[int]:
        """Search for an object

        Args:
            diameter: diameter to search for (optional)
            attribute: attribute to search for (optional)

        Returns:
            indices of found objects
        """
        rtn = []

        for i, _ in enumerate(self._np_shapes.diameter):
            if (diameter is not None and self._np_shapes.diameter[i] != diameter) or \
                    (attribute is not None and
                        self._attributes[i] != attribute):
                continue
            rtn.append(i)
        return rtn

    def csv(self, variable_names: bool = True,
            hash_column: bool = False,
            attribute_column: bool = True) -> str:
        # inherited doc
        size_dict = {"diameter": self._np_shapes.diameter}
        if attribute_column:
            attr = self.attributes
        else:
            attr = None
        if hash_column:
            array_hash = self.hash
        else:
            array_hash = None
        return make_csv(xy=self._np_shapes.xy,
                        size_data_dict=size_dict,
                        attributes=attr, array_hash=array_hash,
                        make_variable_names=variable_names)
