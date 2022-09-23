from typing import Optional

import numpy as np
from numpy.typing import ArrayLike
from .._lib.typing import IntOVector
from .._lib.np_tools import make_vector_fixed_length

from .dot import Dot
from .rectangle import Rectangle


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
