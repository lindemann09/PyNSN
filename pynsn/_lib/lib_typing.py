from __future__ import annotations

from typing import Tuple, Union, Sequence, Optional, List, Any, Iterator
from numpy.typing import ArrayLike, NDArray
from numpy import float_, int_

# Typing
IntOVector = Union[int, Sequence[int], NDArray[int_]]
StrOVector = Union[str, Sequence[str], NDArray[str]]
FloatOVector = Union[float, Sequence[float], NDArray[float_]]
NumPair = Union[Tuple[float, float], List[float], NDArray[float_]]
OptArrayLike = Optional[ArrayLike]
OptInt = Optional[int]
OptStr = Optional[str]
OptFloat = Optional[float]
