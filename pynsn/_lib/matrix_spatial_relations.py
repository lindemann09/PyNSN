"""Numpy optimized function of geometry"""

__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from itertools import combinations

import numpy as np
from numpy.typing import NDArray

from .spatial_relations import CoordinateSpatRel, DotSpatRel, RectangleSpatRel


class _CombinationMatrix(object):
    """Symmetric combination matrix
    helper class"""

    def __init__(self, n_items: int) -> None:
        idx = np.array(list(combinations(range(n_items), r=2)))
        self.idx_a = idx[:, 0]
        self.idx_b = idx[:, 1]
        self.matrix = np.full((n_items, n_items), np.nan)

    @property
    def n_combinations(self) -> float:
        return len(self.idx_a)

    def fill(self, values) -> NDArray:
        """returns combination matrix with values

        """
        self.matrix[self.idx_a, self.idx_b] = values
        self.matrix[self.idx_b, self.idx_a] = values
        return self.matrix


class RectangleMatrixSpatRel(_CombinationMatrix):

    def __init__(self, xy: NDArray, sizes: NDArray):
        super().__init__(n_items=xy.shape[0])
        self._rr = RectangleSpatRel(a_xy=xy[self.idx_a, :],
                                    a_sizes=sizes[self.idx_a, :],
                                    b_xy=xy[self.idx_b, :],
                                    b_sizes=sizes[self.idx_b, :])

    def overlaps(self, minimum_distance: float = 0) -> NDArray:
        """Matrix with overlaps (True/False) between the rectangles"""
        return self.fill(values=self._rr.overlaps(minimum_distance))

    def distances(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self.fill(values=self._rr.distances())


class CoordinateMatrixSpatRel(_CombinationMatrix):

    def __init__(self, xy: NDArray):
        super().__init__(n_items=xy.shape[0])
        self._rr = CoordinateSpatRel(a_xy=xy[self.idx_a, :],
                                     b_xy=xy[self.idx_b, :])

    def distances(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self.fill(values=self._rr.distances())


class DotMatrixSpatRel(_CombinationMatrix):

    def __init__(self, xy: NDArray, diameter: NDArray):
        super().__init__(n_items=xy.shape[0])
        self._rr = DotSpatRel(a_xy=xy[self.idx_a, :],
                              a_diameter=diameter[self.idx_a, :],
                              b_xy=xy[self.idx_b, :],
                              b_diameter=diameter[self.idx_b, :])

    def overlaps(self, minimum_distance: float = 0) -> NDArray:
        """Matrix with overlaps (True/False) between the rectangles"""
        return self.fill(values=self._rr.overlaps(minimum_distance))

    def distances(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self.fill(values=self._rr.distances())
