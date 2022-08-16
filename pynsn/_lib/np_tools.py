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


def as_vector(x: ArrayLike) -> NDArray:
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


def make_array2d(arr: NDArray, n_rows: int) -> NDArray:
    """make 2d array with n equal rows, if array is zero or one dimensional"""
    if arr.ndim < 2 or arr.shape[0] != n_rows:
        try:
            return np.ones((n_rows, 1)) * arr
        except ValueError as err:
            raise ValueError(f"Can not make a numpy array with {n_rows} rows from "
                             f" an array with the shape={arr.shape}.") from err
    else:
        return arr


def as_array2d(arr: ArrayLike) -> NDArray:
    """converts a simple 1D ArrayLike object to a 2D numpy array with only
    one row. It thus insures a multi-dimensional numpy array"""
    rtn = np.asarray(arr)
    if rtn.ndim == 1:
        return rtn.reshape((1, rtn.shape[0]))
    return rtn


def round2(arr: NDArray, decimals: int,
           int_type: type = np.int32) -> NDArray:
    """rounds and changes to int type if decimals == 0"""
    arr = np.round(arr, decimals=decimals)
    if decimals == 0:
        return arr.astype(int_type)
    else:
        return arr


def is_all_equal(vector: ArrayLike)->bool:
    # returns true if all elements are equal
    return len(np.unique(np.asarray(vector))) == 1


def triu_nan(m: NDArray, k:int=0)-> NDArray:
    """helper function
    upper triangular but nan instead of zeros (as in numpy's original function,
    see docu numpy.triu)
    """
    return m + np.tril(np.full(m.shape, np.nan), k=k-1)
