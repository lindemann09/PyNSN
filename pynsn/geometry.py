from copy import deepcopy as _deepcopy

import numpy as _np
from numpy.typing import NDArray as _NDArray

from ._stimulus.nsn_stimulus import NSNStimulus as _NSNStimulus


def distances(shape_array: _NSNStimulus) -> _NDArray:
    """Matrix with the distances between all shapes in the array"""
    return _relation_matrix(shape_array, what=0)


def dwithin(shape_array: _NSNStimulus, distance: float) -> _NDArray:
    """Matrix with the booleans indicates dwithin between all shapes in the array"""
    return _relation_matrix(shape_array, what=1, para=distance)


def overlaps(shape_array: _NSNStimulus, min_distance: float | None = None) -> _NDArray:
    if min_distance is None:
        min_distance = shape_array.min_distance
    else:
        min_distance = 0

    return dwithin(shape_array, distance=min_distance)


def _relation_matrix(arr: _NSNStimulus, what: int, para: float = 0) -> _NDArray:
    """helper function returning the relation between polygons
    0 = distance
    1 = dwithin
    """
    arr = _deepcopy(arr)
    l = arr.n_shapes
    rtn = _np.full((l, l), _np.nan)
    for x in reversed(range(l)):
        shape = arr.shape_pop(x)
        if what == 0:
            y = arr.distances(shape)
        elif what == 1:
            y = arr.dwithin(shape=shape, distance=para)
        else:
            raise RuntimeError("unknown function")

        rtn[x, 0:x] = y

    # make symmetric
    i_lower = _np.triu_indices(l, 1)
    rtn[i_lower] = rtn.T[i_lower]
    return rtn
