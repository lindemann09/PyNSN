__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from itertools import combinations

import numpy as np
from numpy.typing import ArrayLike, NDArray


class CombinationMatrx(object):
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


def as_vector(x):
    """helper function:
    make an numpy vector from any element (list, _arrays, and single data (str, numeric))
    """

    x = np.asarray(x)
    if x.ndim == 1:
        return x
    elif x.ndim == 0:
        # if one element only, make a array with one element
        return x.reshape(1)
    else:
        return x.flatten()


def as_array2d(data: ArrayLike) -> NDArray[np.floating]:
    """converts a simple 1D array to 2D with only row.
    It just insures a multi-dimensional array"""
    rtn = np.asarray(data)
    if rtn.ndim == 1:
        return rtn.reshape((1, rtn.shape[0]))
    return rtn


def round2(array: NDArray, decimals: int,
           int_type: type = np.int32) -> NDArray:
    """rounds and changes to int type if decimals == 0"""
    array = np.round(array, decimals=decimals)
    if decimals == 0:
        return array.astype(int_type)
    else:
        return array


def is_all_equal(vector):
    # returns true if all elements are equal
    return len(np.unique(np.asarray(vector))) == 1


def triu_nan(m, k=0):
    """helper function
    upper triangular but nan instead of zeros (as in numpy's original function,
    see docu numpy.triu)
    """
    return m + np.tril(np.full(m.shape, np.nan), k=k-1)
