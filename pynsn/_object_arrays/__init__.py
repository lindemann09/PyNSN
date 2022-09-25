import typing as _tp

from .dot_array import BaseDotArray, DotArray
from .rectangle_array import BaseRectangleArray, RectangleArray
from ..spatial_relations import _spatial_relations
from ..spatial_relations import _matrix_spatial_relations

ObjectArrayType = _tp.Union[DotArray, RectangleArray]
