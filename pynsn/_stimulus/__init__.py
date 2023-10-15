import typing as _tp
from .shapes import Dot, Rectangle, Picture
from .properties import ArrayProperties, VisProp
from .shape_array import ShapeArray
from .nsn_stimulus import NSNStimulus

ShapeType = _tp.Union[Dot, Rectangle, Picture]
