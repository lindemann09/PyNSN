import typing as _tp
from .shapes import Dot, Rectangle, Picture
from .properties import ArrayProperties
from .shape_array import ShapeArray
from .nsn_stimulus import NSNStimulus

ShapeType = _tp.Union[Dot, Rectangle, Picture]
