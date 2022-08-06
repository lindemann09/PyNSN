import typing as _tp

from .load_array import load_array
from .properties import ArrayProperties
from .rect_array import RectangleArray
from .dot_array import DotArray
from .point_array import PointArray
from .target_area import TargetArea

# TODO not yet used
ObjectArrayType = _tp.Union[DotArray, RectangleArray]
