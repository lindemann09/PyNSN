import typing as _tp

from .dot import Dot
from .rectangle import Rectangle
from .picture import Picture
from .coordinate import Coordinate

ShapeType = _tp.Union[Dot, Rectangle, Picture]
RectangleLike = _tp.Union[Rectangle, Picture]
