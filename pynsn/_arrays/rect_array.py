"""
Rectangle Array
"""
from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from copy import deepcopy
from typing import Any, Iterator, Optional, Sequence, Union

import numpy as np
from numpy.typing import NDArray, ArrayLike

from pynsn._shapes import ShapeType

from .._lib import np_tools
from .._lib.coordinate import Coordinate
from .._lib.spatial_relations import RectangleRectangle
from .._lib.typing import IntOVector
from .._shapes.rectangle import Rectangle
from .abc_object_array import ABCObjectArray
from .tools import make_csv
from .._shapes.rectangle_list import RectangleList

# pylint: disable=W0237:arguments-renamed


class RectangleArray(ABCObjectArray):
    """Rectangle array class

    All object arrays are restricted to a certain target area. The classes are
    optimized for numpy calculations
    """

    def __init__(self,
                 target_area: ShapeType,
                 min_distance_between_objects: Optional[float] = None,
                 min_distance_area_boarder: Optional[float] = None,
                 xy: Optional[ArrayLike] = None,
                 sizes: Optional[ArrayLike] = None,
                 attributes: Optional[ArrayLike] = None) -> None:
        """Rectangular array is restricted to a certain area, it has a target area
        and a minimum gap.

        This properties allows shuffling free position and adapting
        properties.

        """
        self._objs = RectangleList(xy=xy, sizes=sizes)
        super().__init__(object_shapes=self._objs,
                         attributes=attributes,
                         target_area=target_area,
                         min_distance_between_objects=min_distance_between_objects,
                         min_distance_area_boarder=min_distance_area_boarder)

    def add(self, obj: Union[Rectangle, Sequence[Rectangle]]) -> None:
        """append one dot or list of dots"""
        if isinstance(obj, Rectangle):
            obj = [obj]
        else:
            obj = list(obj)

        xy = []
        sizes = []
        attributes = []
        for r in obj:
            xy.append(r.xy)
            sizes.append(r.size)
            attributes.append(r.attribute)

        self._objs.append(xy, sizes)
        self._set_missing_attributes(attributes)
        self.properties.reset()

    @property
    def sizes(self) -> NDArray:
        """Numpy array of the object sizes

        The two dimensional array (shape=[2, n]) represents the width and height of the `n` objects in this array

        Returns:
            the width and height

        """
        return self._objs.sizes

    @property
    def surface_areas(self) -> NDArray:
        # a = w*h
        return self._objs.sizes[:, 0] * self._objs.sizes[:, 1]

    @property
    def perimeter(self) -> NDArray:
        return 2 * (self._objs.sizes[:, 0] + self._objs.sizes[:, 1])

    def mod_round_values(self, decimals: int = 0,
                         int_type: type = np.int16) -> None:
        # inherited doc
        if decimals is None:
            return
        self._objs.xy = np_tools.round2(self._objs.xy, decimals=decimals,
                                        int_type=int_type)
        self._objs.sizes = np_tools.round2(self._objs.sizes, decimals=decimals,
                                           int_type=int_type)

    def to_dict(self) -> dict:
        # inherited doc
        d = super().to_dict()
        d.update({"sizes": self._objs.sizes.tolist()})
        return d

    @staticmethod
    def from_dict(the_dict: dict) -> RectangleArray:
        """read rectangle array from dict"""

        rtn = RectangleArray(target_area=RectangleArray._target_area_from_dict(the_dict),
                             min_distance_between_objects=the_dict["min_distance_between_objects"],
                             min_distance_area_boarder=the_dict["min_distance_area_boarder"])

        rtn._objs.append(xy=the_dict["xy"], sizes=the_dict["sizes"])
        rtn._set_missing_attributes(attributes=the_dict["attributes"])
        return rtn

    def dataframe_dict(self, hash_column: bool = False,
                       attribute_column: bool = True) -> dict:
        # inherited doc
        if hash_column:
            d = {"hash": [self.hash] * len(self._objs.xy)}
        else:
            d = {}
        d.update({"x": self._objs.xy[:, 0].tolist(),
                  "y": self._objs.xy[:, 1].tolist(),
                  "width": self._objs.sizes[:, 0].tolist(),
                  "height": self._objs.sizes[:, 1].tolist()})
        if attribute_column:
            d.update({"attributes": self._attributes.tolist()})
        return d

    def copy(self, indices: Union[int, Sequence[int], None] = None,
             deep_copy: bool = True) -> RectangleArray:
        """returns a (deep) copy of the dot array.

        It allows to copy a subset of dot only.

        """

        if len(self._objs.xy) == 0:
            return RectangleArray(
                target_area=deepcopy(self.target_area),
                min_distance_area_boarder=self.min_distance_area_boarder,
                min_distance_between_objects=self.min_distance_between_objects)
        if indices is None:
            indices = list(range(len(self._objs.xy)))

        if deep_copy:
            return RectangleArray(
                target_area=deepcopy(self.target_area),
                min_distance_between_objects=self.min_distance_between_objects,
                min_distance_area_boarder=self.min_distance_area_boarder,
                xy=self._objs.xy[indices, :].copy(),
                sizes=self._objs.sizes[indices].copy(),
                attributes=self._attributes[indices].copy())
        else:
            return RectangleArray(
                target_area=deepcopy(self.target_area),
                min_distance_between_objects=self.min_distance_between_objects,
                min_distance_area_boarder=self.min_distance_area_boarder,
                xy=self._objs.xy[indices, :],
                sizes=self._objs.sizes[indices],
                attributes=self._attributes[indices])

    def get_distances(self, rect: Rectangle) -> NDArray:
        """Euclidean distances toward a single rectangle

        Negative numbers indicate an overlap of the objects

        Args:
            rect: object towards the distances will be calculated

        Returns: numpy array of distances
        """

        assert isinstance(rect, Rectangle)
        if len(self._objs.xy) == 0:
            return np.array([])
        else:
            rel = RectangleRectangle(a_xy=self._objs.xy,
                                     a_sizes=self._objs.sizes,
                                     b_xy=rect.xy,
                                     b_sizes=rect.size)
            return rel.xy_distances()

    def iter_objects(self, indices: Optional[IntOVector] = None) -> Iterator[Rectangle]:
        """iterate over all or a part of the objects

        Parameters
        ----------
        indices

        Notes
        -----
        To iterate all object you might all use the class iterator __iter__:
        >>> for obj in my_array:
        >>>    print(obj)
        """

        if isinstance(indices, (int, np.integer)):
            yield Rectangle(xy=self._objs.xy[indices, :],
                            size=self._objs.sizes[indices],
                            attribute=self._attributes[indices])
        else:
            if indices is None:
                data = zip(self._objs.xy,
                           self._objs.sizes, self._attributes)
            else:
                data = zip(self._objs.xy[indices, :],  # type: ignore
                           self._objs.sizes[indices, :],  # type: ignore
                           self._attributes[indices])

            for xy, s, att in data:
                rtn = Rectangle(xy=xy, size=s, attribute=att)
                yield rtn

    def find_objects(self, size: Union[Sequence[float], NDArray[np.float_], None] = None,
                     attribute: Optional[Any] = None,
                     corner: Optional[Coordinate] = None) -> Sequence[int]:
        """returns indices of found objects

        2D-tuple
        """
        rtn = []
        for i in range(len(self._objs.sizes)):
            if (size is not None and
                self._objs.sizes[i, 0] != size[0] and self._objs.sizes[i, 1] != size[1]) or \
                    (attribute is not None and self._attributes[i] != attribute):
                continue
            rtn.append(i)

        if corner is None:
            return rtn
        elif isinstance(corner, Coordinate):
            new_rtn = []
            for i, rect in zip(rtn, self.iter_objects(indices=rtn)):
                if corner in list(rect.iter_corners()):
                    new_rtn.append(i)
            return new_rtn
        else:
            raise TypeError("corner has to be a Coordinate and not {}.".format(
                type(corner)))

    def csv(self, variable_names: bool = True,
            hash_column: bool = False,
            attribute_column: bool = True) -> str:
        # inherited doc
        size_dict = {
            "width": self._objs.sizes[:, 0], "height": self._objs.sizes[:, 1]}
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
