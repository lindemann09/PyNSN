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

    @property
    def n_combinations(self) -> float:
        return len(self.idx_a)

    def get_matrix(self, values) -> NDArray:
        """returns matrix with [idx_a, idx_b, values, ...]

        """
        idx = np.array([self.idx_a, self.idx_b]).T
        if values.ndim == 1:
            values = values.reshape((len(values), 1))
        return np.append(idx, values, axis=1)

    def get_squared_matrix(self, values, both_triangle=True) -> NDArray:
        """returns combination matrix with values

        """
        rtn = np.full((len(self.idx_a), len(self.idx_a)), np.nan)
        rtn[self.idx_a, self.idx_b] = values
        if both_triangle:
            rtn[self.idx_b, self.idx_a] = values
        return rtn


class RectangleMatrixSpatRel(_CombinationMatrix):

    def __init__(self, xy: NDArray, sizes: NDArray):
        super().__init__(n_items=xy.shape[0])
        self._rr = RectangleSpatRel(a_xy=xy[self.idx_a, :],
                                    a_sizes=sizes[self.idx_a, :],
                                    b_xy=xy[self.idx_b, :],
                                    b_sizes=sizes[self.idx_b, :])

    def overlaps(self, minimum_distance: float = 0) -> NDArray:
        """Matrix with overlaps (True/False) between the rectangles"""
        return self.get_matrix(values=self._rr.overlaps(minimum_distance))

    def distances(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self.get_matrix(values=self._rr.distances())

    def spatial_relations(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self.get_matrix(values=self._rr.spatial_relations())

    def required_displacements(self, minimum_distance: float = 0) -> NDArray:
        """Return matrix with distance between the rectangles"""
        rtn = self.get_matrix(
            values=self._rr.required_displacements(minimum_distance))
        i = np.flatnonzero(rtn[:, 2] != 0)
        return rtn[i, :]


class CoordinateMatrixSpatRel(_CombinationMatrix):

    def __init__(self, xy: NDArray):
        super().__init__(n_items=xy.shape[0])
        self._rr = CoordinateSpatRel(a_xy=xy[self.idx_a, :],
                                     b_xy=xy[self.idx_b, :])

    def distances(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self.get_matrix(values=self._rr.distances())

    def spatial_relations(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self.get_matrix(values=self._rr.spatial_relations())


class DotMatrixSpatRel(_CombinationMatrix):

    def __init__(self, xy: NDArray, diameter: NDArray):
        super().__init__(n_items=xy.shape[0])
        self._rr = DotSpatRel(a_xy=xy[self.idx_a, :],
                              a_diameter=diameter[self.idx_a, :],
                              b_xy=xy[self.idx_b, :],
                              b_diameter=diameter[self.idx_b, :])

    def overlaps(self, minimum_distance: float = 0) -> NDArray:
        """Matrix with overlaps (True/False) between the rectangles"""
        return self.get_matrix(values=self._rr.overlaps(minimum_distance))

    def distances(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self.get_matrix(values=self._rr.distances())

    def spatial_relations(self) -> NDArray:
        """Return matrix with distance between the rectangles"""
        return self.get_matrix(values=self._rr.spatial_relations())

    def required_displacements(self, minimum_distance: float = 0) -> NDArray:
        """Return matrix with distance between the rectangles"""
        rtn = self.get_matrix(
            values=self._rr.required_displacements(minimum_distance))
        i = np.flatnonzero(rtn[:, 2] != 0)
        return rtn[i, :]
