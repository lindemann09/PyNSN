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
from .._lib.typing import IntOVector
from .._shapes import Dot,  ShapeType
from .._shapes.dot_list import DotList
from .abc_object_array import ABCObjectArray
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

        self._objs = DotList(xy, diameter, attributes)  # for typing
        super().__init__(object_list=self._objs,
                         target_area=target_area,
                         min_distance_between_objects=min_distance_between_objects,
                         min_distance_area_boarder=min_distance_area_boarder)

    def add(self, shapes: Union[Dot, Sequence[Dot]]) -> None:
        """append one dot or list of dots"""
        self._objs.add_objects(shapes)

    @property
    def diameter(self) -> NDArray:
        return self._objs.diameter

    def to_dict(self) -> dict:
        d = super().to_dict()
        d.update({"diameter": self._objs.diameter.tolist()})
        return d

    def dataframe_dict(self, hash_column: bool = False,
                       attribute_column: bool = True) -> dict:
        # inherited doc
        if hash_column:
            d = {"hash": [self.hash] * len(self._objs.xy)}
        else:
            d = {}
        d.update({"x": self._objs.xy[:, 0].tolist(),
                  "y": self._objs.xy[:, 1].tolist(),
                  "diameter": self._objs.diameter.tolist()})
        if attribute_column:
            d.update({"attributes": self._objs.attributes.tolist()})
        return d

    @staticmethod
    def from_dict(the_dict: Dict[str, Any]) -> DotArray:
        """read Dot collection from dict"""
        return DotArray(xy=the_dict["xy"],
                        diameter=the_dict["diameter"],
                        attributes=the_dict["attributes"],
                        target_area=DotArray._target_area_from_dict(the_dict),
                        min_distance_between_objects=the_dict["min_distance_between_objects"],
                        min_distance_area_boarder=the_dict["min_distance_area_boarder"])

    def copy(self, indices: Union[int, Sequence[int], NDArray[np.int_], None] = None,
             deep_copy: bool = True) -> DotArray:
        """A (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        Args:
            indices: Arrays with the indices to be be copied. If None, all
                object will be copies

        Returns:
            a copy of the array
        """

        if len(self._objs.xy) == 0:
            return DotArray(target_area=deepcopy(self.target_area),
                            min_distance_between_objects=self.min_distance_between_objects,
                            min_distance_area_boarder=self.min_distance_area_boarder)

        if indices is None:
            indices = np.arange(len(self._objs.xy.shape))

        if deep_copy:
            return DotArray(target_area=deepcopy(self.target_area),
                            min_distance_between_objects=self.min_distance_between_objects,
                            min_distance_area_boarder=self.min_distance_area_boarder,
                            xy=self._objs.xy[indices, :].copy(),
                            diameter=self._objs.diameter[indices].copy(),
                            attributes=self._objs.attributes[indices].copy())
        else:
            return DotArray(target_area=deepcopy(self.target_area),
                            min_distance_between_objects=self.min_distance_between_objects,
                            min_distance_area_boarder=self.min_distance_area_boarder,
                            xy=self._objs.xy[indices, :],
                            diameter=self._objs.diameter[indices],
                            attributes=self._objs.attributes[indices])

    def get_distances(self, dot: Dot) -> NDArray:
        """Euclidean distances toward a single dot

        Negative numbers indicate an overlap of the objects

        Args:
            dot: object towards the distances will be calculated

        Returns: numpy array of distances
        """
        assert isinstance(dot, Dot)
        if len(self._objs.xy) == 0:
            return np.array([])
        else:
            rel = DotDot(a_xy=self._objs.xy,
                         a_diameter=self._objs.diameter,
                         b_xy=dot.xy,
                         b_diameter=np.array([dot.diameter]))
            return rel.distances()

    def iter_objects(self, indices: Optional[IntOVector] = None) -> Iterator[Dot]:
        # inherited doc
        return self._objs.iter_objects(indices)

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

        for i, _ in enumerate(self._objs.diameter):
            if (diameter is not None and self._objs.diameter[i] != diameter) or \
                    (attribute is not None and
                        self._objs.attributes[i] != attribute):
                continue
            rtn.append(i)
        return rtn

    def csv(self, variable_names: bool = True,
            hash_column: bool = False,
            attribute_column: bool = True) -> str:
        # inherited doc
        size_dict = {"diameter": self._objs.diameter}
        if attribute_column:
            attr = self.attributes
        else:
            attr = None
        if hash_column:
            array_hash = self.hash
        else:
            array_hash = None
        return make_csv(xy=self._objs.xy,
                        size_data_dict=size_dict,
                        attributes=attr, array_hash=array_hash,
                        make_variable_names=variable_names)
