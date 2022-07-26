from __future__ import annotations

from typing import Tuple, Union, Sequence, Optional, List, Any, Iterator, Iterable
from numpy.typing import NDArray, ArrayLike
from numpy import float_, int_

# Typing
IntOVector = Union[int, Sequence[int]]
StrOVector = Union[str, Sequence[str]]
FloatOVector = Union[float, Sequence[float], NDArray[float_]]
NumPair = Union[Tuple[float, float], Sequence[float], NDArray[float_]]
NumSeq = Union[Sequence[float_], NDArray[float_]]
OptArrayLike = Optional[ArrayLike]
OptInt = Optional[int]
OptStr = Optional[str]
OptFloat = Optional[float]

# FIXME when to use ArrayLike and when NDArray
