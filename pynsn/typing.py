from typing import Optional, Union, List, Sequence

from numpy import float_, int_, str_
from numpy.typing import ArrayLike, NDArray

OptArrayLike = Optional[ArrayLike]
OptInt = Optional[int]
OptStr = Optional[str]
OptFloat = Optional[float]

IntOVector = Union[int, List[int], Sequence[int_], NDArray[int_]]
StrOVector = Union[str, List[str], Sequence[str_], NDArray[str_]]
FloatOVector = Union[float, List[float], Sequence[float_], NDArray[float_]]
