import typing as _tp

from .load_array import load_array
from .properties import ArrayProperties
from .rect_array import RectangleArray
from .dot_array import DotArray
from .target_area import TargetArea

ObjectArrayType = _tp.Union[DotArray, RectangleArray]
