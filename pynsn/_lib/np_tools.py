__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Iterable, List

import numpy as np
from numpy.typing import ArrayLike, NDArray


def as_vector(x: ArrayLike) -> NDArray:
    """helper function:
    make an numpy vector from any element (list, _arrays, and single data (str, numeric))
    """

    return np.atleast_1d(x).flatten()


def as_array2d_nrow(arr: ArrayLike, n_rows: int) -> NDArray:
    """make 2d array with n equal rows, if array is zero or one dimensional"""
    arr = np.atleast_2d(arr)
    if arr.shape[0] != n_rows:
        try:
            return np.ones((n_rows, 1)) * arr
        except ValueError as err:
            raise ValueError(f"Can not make a numpy array with {n_rows} rows from "
                             f" an array with the shape={arr.shape}.") from err
    else:
        return arr


def all_as_array2d_equal_rows(list_of_arrays: Iterable[ArrayLike]) -> List[NDArray]:
    """make all arrays 2D. If one array is 2d with n rows all zero or one
    dimensional arrays will be converted to arrays with n identical rows"""
    n_rows = 1
    array_list = []
    for arr in list_of_arrays:
        array_list.append(np.asarray(arr))
        if array_list[-1].ndim > 1 and array_list[-1].shape[0] > n_rows:
            n_rows = array_list[-1].shape[0]

    for i, arr in enumerate(array_list):
        array_list[i] = as_array2d_nrow(arr, n_rows=n_rows)

    return array_list


def round2(arr: NDArray, decimals: int,
           int_type: type = np.int32) -> NDArray:
    """rounds and changes to int type if decimals == 0"""
    arr = np.round(arr, decimals=decimals)
    if decimals == 0:
        return arr.astype(int_type)
    else:
        return arr


def is_all_equal(vector: ArrayLike) -> np.bool_:
    # returns true if all elements are equal
    a = np.asarray(vector).flatten()
    return np.all(np.equal(a, a[0]))


def triu_nan(m: NDArray, k: int = 0) -> NDArray:
    """helper function
    upper triangular but nan instead of zeros (as in numpy's original function,
    see docu numpy.triu)
    """
    return m + np.tril(np.full(m.shape, np.nan), k=k-1)
