import numpy as _np
import typing as _tp

NumVector = Sequence[float]
NumArray = Sequence[NumVector]
NumVectorOrArray = Union[NumArray, NumVector]
Attribute = Union[str, int, float, None]
AttributeOrAttributeVector = Union[Attribute, Sequence(Attribute)]

NPArray = _np.ndarray