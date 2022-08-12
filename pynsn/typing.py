from typing import Optional, Union, List

from numpy import float_, int_, str_
from numpy.typing import ArrayLike, NDArray

OptArrayLike = Optional[ArrayLike]
OptInt = Optional[int]
OptStr = Optional[str]
OptFloat = Optional[float]

IntOVector = Union[int, List[int], NDArray[int_]]
StrOVector = Union[str, List[str], NDArray[str_]]
FloatOVector = Union[float, List[float], NDArray[float_]]
