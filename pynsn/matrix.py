
from copy import deepcopy
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from ._stimulus.nsn_stimulus import NSNStimulus
from ._stimulus.shape_array import ShapeArray


def distances(shape_array:ShapeArray) -> NDArray:
    """Matrix with the distances between the shapes in the array"""
    return _relation_matrix(shape_array, what=0)

def dwithin(shape_array:ShapeArray, distance: float) -> NDArray:
    """Matrix with the booleans indication dwith between the shapes in the array"""
    return _relation_matrix(shape_array, what=1, para=distance)

def overlaps(shape_array:ShapeArray, min_distance: Optional[float] = None) -> NDArray:
    if min_distance is None:
        if isinstance(shape_array, NSNStimulus):
            min_distance = shape_array.min_distance
        else:
            min_distance = 0

    return dwithin(shape_array, distance=min_distance)

def _relation_matrix(arr:ShapeArray, what: int, para: float = 0) -> NDArray:
    """helper function returning the relation between polygons
    0 = distance
    1 = dwithin
    """
    arr = deepcopy(arr)
    l = arr.n_objects
    rtn = np.full((l, l), np.nan)
    for x in reversed(range(l)):
        shape = arr.pop(x)
        if what == 0:
            y = arr.distances(shape)
        elif what == 1:
            y = arr.dwithin(shape=shape, distance=para)
        else:
            raise RuntimeError("unknown function")

        rtn[x, 0:x] = y

    # make symmetric
    i_lower = np.triu_indices(l, 1)
    rtn[i_lower] = rtn.T[i_lower]
    return rtn