from __future__ import annotations

from typing import Optional, Sequence, Union

from numpy import float_
from numpy.typing import ArrayLike, NDArray

IntOVector = Union[int, Sequence[int]]
StrOVector = Union[str, Sequence[str]]
FloatOVector = Union[float, Sequence[float], NDArray[float_]]
NumArray = NDArray[float_]
OptArrayLike = Optional[ArrayLike]
OptInt = Optional[int]
OptStr = Optional[str]
OptFloat = Optional[float]
