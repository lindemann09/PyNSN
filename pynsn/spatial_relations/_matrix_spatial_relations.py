"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from itertools import combinations
from typing import Union

import numpy as np
from numpy.typing import NDArray

from .._object_arrays import BaseDotArray, BaseRectangleArray
from ._spatial_relations import DotDot, RectangleRectangle


class SpatRelMatrix():

    def __init__(self, object_array: Union[BaseDotArray, BaseRectangleArray]):

        if not isinstance(object_array, (BaseRectangleArray, BaseDotArray)):
            raise TypeError(f"Spatial relations matrix requires for BaseRectangleArray "
                            f"or BaseDotArray, but not {type(object_array)}.")

        self.n_objects = len(object_array.xy)
        self._ix = np.array(list(combinations(range(self.n_objects), r=2)))

        if isinstance(object_array, BaseRectangleArray):

            rect_a = BaseRectangleArray(xy=object_array.xy[self._ix[:, 0], :],
                                        sizes=object_array.sizes[self._ix[:, 0], :])
            rect_b = BaseRectangleArray(xy=object_array.xy[self._ix[:, 1], :],
                                        sizes=object_array.sizes[self._ix[:, 1], :])
            self._rr = RectangleRectangle(a_rectangles=rect_a,
                                          b_rectangles=rect_b)
        elif isinstance(object_array, BaseDotArray):
            dots_a = BaseDotArray(xy=object_array.xy[self._ix[0], :],
                                  diameter=object_array.diameter[self._ix[0], :])
            dots_b = BaseDotArray(xy=object_array.xy[self._ix[1], :],
                                  diameter=object_array.diameter[self._ix[1], :])
            self._rr = DotDot(a_dots=dots_a, b_dots=dots_b)

    def _matrix(self, values: NDArray) -> NDArray:
        """returns matrix with values"""
        rtn = np.full(shape=(self.n_objects, self.n_objects),
                      fill_value=np.nan)
        rtn[self._ix[:, 0], self._ix[:, 1]] = values
        return rtn

    def is_inside(self, minimum_gap: float = 0) -> NDArray:
        """Return matrix with distance between the objects"""
        return self._matrix(values=self._rr.is_inside(minimum_gap))

    def overlaps(self, minimum_gap: float = 0) -> NDArray:
        """Return matrix with distance between the objects"""
        return self._matrix(values=self._rr.overlaps(minimum_gap))

    def distances(self) -> NDArray:
        """Return matrix with distance between the objects"""
        return self._matrix(values=self._rr.distances)

    def distances_radial(self) -> NDArray:
        return self._matrix(values=self._rr.distances_radial)
