from typing import Optional, Sequence, Union, Iterator

import numpy as np
from numpy.typing import ArrayLike, NDArray
from .._lib.typing import IntOVector
from .._lib.np_tools import make_vector_fixed_length, round2

from .dot import Dot


class BaseDotList:
    """Basic class comprising numpy lists of coordinates and diameter for dots"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 diameter: Optional[ArrayLike] = None) -> None:
        self.xy = np.empty((0, 2))
        self.diameter = np.array([])
        if xy is not None and diameter is not None:
            self.xy = np.atleast_2d(xy)
            self.diameter = np.atleast_1d(diameter)
            self.check()

    def check(self):
        """raises value error if badly shaped data"""
        if self.xy.shape[0] != len(self.diameter):
            raise ValueError("Badly shaped data: " +
                             "xy has not the same length as item_diameter")


class DotList(BaseDotList):
    """Class of numpy representation of dot"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 diameter: Optional[ArrayLike] = None,
                 attributes: Optional[ArrayLike] = None) -> None:
        super().__init__()
        self.attributes = np.array([])
        if xy is not None and diameter is not None:
            self.append(xy, diameter, attributes)

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

    def append(self,
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

    def add_objects(self, shapes: Union[Dot, Sequence[Dot]]) -> None:
        """append one dot or list of dots"""
        if isinstance(shapes, Dot):
            shapes = [shapes]
        else:
            shapes = list(shapes)

        for obj in shapes:
            self.xy = np.append(self.xy, np.atleast_2d(obj.xy), axis=0)
            self.diameter = np.append(self.diameter, obj.diameter, axis=0)
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
        """Round all values

        Args:
            decimals: number of decimal places
            int_type: numpy int type (default=numpy.int16)
        """
        if decimals is None:
            return
        self.xy = round2(self.xy, decimals=decimals, int_type=int_type)
        self.diameter = round2(self.diameter, decimals=decimals,
                               int_type=int_type)

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
