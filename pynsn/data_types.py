import numpy as _np
import typing as _tp

NumVector = _tp.Sequence[float]
NumArray = _tp.Sequence[NumVector]
NumVectorOrArray = _tp.Union[NumArray, NumVector]
Attribute = _tp.Union[str, int, float, None]
AttributeOrAttributeVector = _tp.Union[Attribute, _tp.Sequence(Attribute)]

NPArray = _np.ndarray