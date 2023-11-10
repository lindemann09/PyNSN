from typing import Sequence, Tuple, Union

from numpy import int_
from numpy.typing import NDArray

IntOVector = Union[int, Sequence[int], NDArray[int_]]

Coord2D = Tuple[float, float]
Coord2DLike = Union[Coord2D, Sequence[float], NDArray]


class NoSolutionError(Exception):
    pass
