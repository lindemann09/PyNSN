__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Optional

import numpy as np
from numpy.typing import ArrayLike

from .._lib.typing import IntOVector
from .._lib.np_tools import make_vector_fixed_length

from .dot import Dot
from .rectangle import Rectangle


class BaseRectangleList:
    """Basic class comprising numpy lists of coordinates and sizes for rectangles"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 sizes: Optional[ArrayLike] = None) -> None:
        self.xy = np.empty((0, 2))
        self.sizes = np.empty((0, 2))

        if xy is not None and sizes is not None:
            self.xy = np.atleast_2d(xy)
            self.sizes = np.atleast_2d(sizes)
            self.check()

    def check(self):
        """raises value error if badly shaped data"""
        if self.xy.shape != self.sizes.shape:
            raise ValueError("Badly shaped data: " +
                             "xy has not the same shape as sizes array")


class RectangleList(BaseRectangleList):
    """Class of numpy representation of rectangles"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 sizes: Optional[ArrayLike] = None,
                 attributes: Optional[ArrayLike] = None) -> None:
        super().__init__()
        self.attributes = np.array([])
        if xy is not None and sizes is not None:
            self.append(xy, sizes, attributes)

    def check(self) -> None:
        """raises value error if badly shaped data"""
        super().check()
        if len(self.attributes) != self.xy.shape[0]:
            raise ValueError("Badly shaped data: attributes have not " +
                             "the correct length.")

    def append(self,
               xy: ArrayLike,
               sizes: ArrayLike,
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
        self.sizes = np.append(self.sizes, np.atleast_2d(sizes), axis=0)
        self.attributes = np.append(self.attributes, attributes)
        self.check()

    def delete(self, index: IntOVector) -> None:
        self.xy = np.delete(self.xy, index, axis=0)
        self.sizes = np.delete(self.xy, index, axis=0)
        self.attributes = np.delete(self.attributes, index)

    def clear(self):
        self.xy = np.empty((0, 2))
        self.sizes = np.empty((0, 2))
        self.attributes = np.array([])
