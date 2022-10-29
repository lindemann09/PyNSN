from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional, Sequence, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._lib.np_tools import make_vector_fixed_length, round2
from .._lib.typing import IntOVector
from .._shapes.dot import Dot
from .abc_object_aray import ABCObjectArray, hash_array, make_csv


class BaseDotArray:
    """Basic class comprising numpy lists of coordinates and diameter for dots"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 diameter: Optional[ArrayLike] = None) -> None:
        self.xy = np.empty((0, 2))
        self.diameter = np.array([])
        if xy is not None and diameter is not None:
            self.xy = np.atleast_2d(xy)
            self.diameter = np.atleast_1d(diameter).astype(float)
            self.check()

    def check(self) -> None:
        if self.xy.shape[0] != len(self.diameter):
            raise ValueError("Badly shaped data: " +
                             "xy has not the same length as item_diameter")

    def n_objects(self) -> int:
        return len(self.xy)


class DotArray(BaseDotArray, ABCObjectArray):
    """Class of numpy representation of dot"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 diameter: Optional[ArrayLike] = None,
                 attributes: Optional[ArrayLike] = None) -> None:
        super().__init__()
        self.attributes = np.array([])
        if xy is not None and diameter is not None:
            self.np_append(xy, diameter, attributes)

    def check(self) -> None:
        """raises value error if badly shaped data"""
        super().check()
        if len(self.attributes) != self.xy.shape[0]:
            raise ValueError("Badly shaped data: attributes have not " +
                             "the correct length.")

    def set_attributes(self, attributes: Optional[ArrayLike]) -> None:
        """Set all attributes

        Args:
            attributes: single attribute or list of attributes
        """
        try:
            self.attributes = make_vector_fixed_length(
                attributes, length=self.xy.shape[0])
        except ValueError as err:
            raise ValueError("Length of attribute list does not match the " +
                             "size of the array.") from err

    def hash(self) -> str:
        return hash_array(xy=self.xy, perimeter=self.perimeter,
                          attributes=self.attributes)

    def np_append(self,
                  xy: ArrayLike,
                  diameter: ArrayLike,
                  attributes: Optional[ArrayLike] = None) -> None:
        """if attributes comprises only one element all new objects will get
        this attribute """
        xy = np.atleast_2d(xy)
        try:
            attributes = make_vector_fixed_length(attributes,
                                                  length=xy.shape[0])
        except ValueError as err:
            raise ValueError("Badly shaped data: attributes have not " +
                             "the correct length.") from err

        self.xy = np.append(self.xy, xy, axis=0)
        self.diameter = np.append(self.diameter,
                                  np.atleast_1d(diameter), axis=0)
        self.attributes = np.append(self.attributes, attributes)
        self.check()

    def delete(self, index: IntOVector) -> None:
        self.xy = np.delete(self.xy, index, axis=0)
        self.diameter = np.delete(self.diameter, index)
        self.attributes = np.delete(self.attributes, index)

    def clear(self):
        self.xy = np.empty((0, 2))
        self.diameter = np.array([])
        self.attributes = np.array([])

    def add(self, shapes: Union[Dot, Sequence[Dot], List[Dot]]) -> None:
        """append one dot or list of dots"""
        if isinstance(shapes, Dot):
            shapes = [shapes]
        else:
            shapes = list(shapes)

        for obj in shapes:
            self.xy = np.append(self.xy, np.atleast_2d(obj.xy), axis=0)
            self.diameter = np.append(self.diameter, obj.diameter)
            self.attributes = np.append(self.attributes, obj.attribute)

    @property
    def surface_areas(self) -> NDArray:
        """TODO

        """
        # a = pi r**2 = pi d**2 / 4
        return np.pi * (self.diameter ** 2) / 4.0

    @property
    def perimeter(self) -> NDArray:
        """Perimeter for each dot
        """
        return np.pi * self.diameter

    def round_values(self, decimals: int = 0,
                     int_type: type = np.int32) -> None:
        if decimals is None:
            return
        self.xy = round2(self.xy, decimals=decimals, int_type=int_type)
        self.diameter = round2(self.diameter, decimals=decimals,
                               int_type=int_type)

    def iter(self, indices: Optional[IntOVector] = None) -> Iterator[Dot]:
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
            yield Dot(xy=self.xy[indices, :],
                      diameter=self.diameter[indices],
                      attribute=self.attributes[indices])
        else:
            if indices is None:
                data = zip(self.xy,
                           self.diameter,
                           self.attributes)
            else:
                data = zip(self.xy[indices, :],
                           self.diameter[indices],
                           self.attributes[indices])
            for xy, dia, att in data:
                yield Dot(xy=xy, diameter=dia, attribute=att)

    def find(self, diameter: Optional[float] = None,
             attribute: Any = None) -> List[int]:
        """returns indices of found objects

        2D-tuple
        """

        rtn = []

        for i, _ in enumerate(self.diameter):
            if (diameter is not None and self.diameter[i] != diameter) or \
                    (attribute is not None and
                        self.attributes[i] != attribute):
                continue
            rtn.append(i)
        return rtn

    def to_dict(self) -> dict:
        d = {"xy": self.xy.tolist()}
        d.update({"diameter": self.diameter.tolist()})
        if len(self.attributes) > 0:
            if np.all(self.attributes == self.attributes[0]):
                # all equal
                d.update({"attributes": self.attributes[0]})
            else:
                d.update({"attributes": self.attributes.tolist()})
        return d

    @classmethod
    def from_dict(cls, the_dict: Dict[str, Any]) -> DotArray:
        """read Dot collection from dict"""
        return cls(xy=the_dict["xy"],
                   diameter=the_dict["diameter"],
                   attributes=the_dict["attributes"])

    def copy(self,
             indices: Union[int, Sequence[int], NDArray[np.int_], None] = None,
             deep_copy: bool = True) -> DotArray:
        # inherited docs
        if indices is None or len(self.xy) == 0:
            xy = self.xy
            diameter = self.diameter
            attributes = self.attributes
        else:
            xy = self.xy[indices, :]
            diameter = self.diameter[indices]
            attributes = self.attributes[indices]

        if deep_copy:
            return DotArray(xy=xy.copy(),
                            diameter=diameter.copy(),
                            attributes=attributes.copy())
        else:
            return DotArray(xy=xy, diameter=diameter, attributes=attributes)

    def dataframe_dict(self, hash_column: bool = False,
                       attribute_column: bool = True) -> dict:
        # inherited docs

        if hash_column:
            d = {"hash": [self.hash()] * len(self.xy)}
        else:
            d = {}
        d.update({"x": self.xy[:, 0].tolist(),
                  "y": self.xy[:, 1].tolist(),
                  "diameter": self.diameter.tolist()})
        if attribute_column:
            d.update({"attributes": self.attributes.tolist()})
        return d

    def csv(self,
            variable_names: bool = True,
            hash_column: bool = False,
            attribute_column: bool = True) -> str:
        # inherited docs

        if attribute_column:
            attr = self.attributes
        else:
            attr = None
        if hash_column:
            array_hash = self.hash()
        else:
            array_hash = None

        return make_csv(xy=self.xy,
                        size_data_dict={"diameter": self.diameter},
                        attributes=attr,
                        array_hash=array_hash,
                        make_variable_names=variable_names)

    @property
    def center_of_mass(self) -> NDArray:
        """center of mass of all objects"""
        weighted_sum = np.sum(
            self.xy * np.atleast_2d(self.perimeter).T, axis=0)
        return weighted_sum / np.sum(self.perimeter)
