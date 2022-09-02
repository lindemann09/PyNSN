import typing as _tp

from .dot import Dot
from .rectangle import Rectangle
from .picture import Picture

ShapeType = _tp.Union[Dot, Rectangle, Picture]
