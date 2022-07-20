from typing import Tuple, Union, Sequence, Optional, List, Any, Iterable
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
ObjectArray = Any # TODO can't be used because of circular import Python 3.7+ see https://stackoverflow.com/questions/42845972/typed-python-using-the-classes-own-type-inside-class-definition
