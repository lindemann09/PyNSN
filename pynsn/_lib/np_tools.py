__author__ = 'Oliver Lindemann <lindemann@cognitive-psychology.eu>'

from typing import Iterable, List, Optional, Tuple

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


def triu_nan(arr: NDArray, k: int = 0) -> NDArray:
    """helper function
    upper triangular but nan instead of zeros (as in numpy's original function,
    see docu numpy.triu)
    """
    return arr + np.tril(np.full(arr.shape, np.nan), k=k-1)


def find_rowwise_one_true(boolean_array2d: NDArray) -> Tuple[NDArray, NDArray]:
    """returns list of the indices (row, column) with one true cell in each row.
    That is, if more than one cell in a row contains True, the function select
    randomly one

    returns idx_row, idx_column
    """
    # (n=rect_id, 2) array of index  (one row -> one cell)
    idx_value = np.vstack(np.nonzero(boolean_array2d)).T
    # Problem: one row could have the values multiple times
    #   --> choose randomly one
    np.random.shuffle(idx_value)  # shuffle rows
    # return only first unique id (ui-index) (Note: set is sorted)
    _, ui_idx = np.unique(idx_value[:, 0], return_index=True)
    return idx_value[ui_idx, 0], idx_value[ui_idx, 1]


def make_vector_fixed_length(values: Optional[ArrayLike], length: int):
    """Helper function: (a) makes vector and scales to n elements, if attributes
    comprises only one element and (b) converts None to nparray([np.nan])
    """

    if values is None:
        values = np.atleast_1d([np.nan])
    else:
        values = np.atleast_1d(values)

    if len(values) == 1 and length > 1:
        values = np.repeat(values, length)
    elif len(values) != length:
        raise ValueError("Values have to be a length of 1 "
                         f"or n_elements ({length})")

    return values


def abs_maximum(arr, axis=None):
    """returns the values with the maximum absolute values

    Note: function avoids call of numpy.abs() for efficiency reasons
    https://stackoverflow.com/questions/17794266/how-to-get-the-highest-element-in-absolute-value-in-a-numpy-matrix
    """
    amax = np.max(arr, axis=axis)
    amin = np.min(arr, axis=axis)
    return np.where(-amin > amax, amin, amax)


def salted_rows(array: NDArray, salt=1e-30) -> NDArray:
    """if array (n, 2) has to identical cells in a row, add salt to one
    randomly chosen cell in that row
    """
    rnd_column = np.random.randint(2, size=len(array))
    array[:, rnd_column] = array[:, rnd_column] + \
        (array[:, 0] == array[:, 1]) * salt
    return array
