from __future__ import annotations

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import List, Optional, Sequence, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .dot import Dot
from .rectangle import Rectangle

IntOVector = Union[int, List[int], Sequence[np.int_], NDArray[np.int_]]


class NPDots(object):
    """Numpy Dots Representation"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 diameter: Optional[ArrayLike] = None) -> None:
        self.xy = np.empty((0, 2))
        self.diameter = np.array([])
        if xy is not None and diameter is not None:
            self.append(xy, diameter)

    def append(self, xy: ArrayLike, diameter: ArrayLike) -> None:
        self.xy = np.append(self.xy, np.atleast_2d(xy), axis=0)
        self.diameter = np.append(self.diameter,
                                  np.atleast_1d(diameter))

        if self.xy.shape[0] != len(self.diameter):
            raise ValueError("Bad shaped data: " +
                             "xy has not the same length as item_diameter")

    def delete(self, index: IntOVector) -> None:
        self.xy = np.delete(self.xy, index, axis=0)
        self.diameter = np.delete(self.diameter, index)

    def clear(self):
        self.xy = np.empty((0, 2))
        self.diameter = np.array([])


class NPRectangles(object):
    """Numpy Rectangle Representation"""

    def __init__(self, xy: Optional[ArrayLike] = None,
                 sizes: Optional[ArrayLike] = None) -> None:
        self.xy = np.empty((0, 2))
        self.sizes = np.empty((0, 2))
        if xy is not None and sizes is not None:
            self.append(xy, sizes)

    def append(self, xy: ArrayLike, sizes: ArrayLike) -> None:
        self.xy = np.append(self.xy, np.atleast_2d(xy), axis=0)
        self.sizes = np.append(self.sizes, np.atleast_2d(sizes), axis=0)

        if self.xy.shape != self.sizes.shape:
            raise ValueError("Bad shaped data: " +
                             "xy has not the same length as sizes array")

    def delete(self, index: IntOVector) -> None:
        self.xy = np.delete(self.xy, index, axis=0)
        self.sizes = np.delete(self.xy, index, axis=0)

    def clear(self):
        self.xy = np.empty((0, 2))
        self.sizes = np.empty((0, 2))
