import typing as _tp

from .dot import Dot
from .rectangle import Rectangle
from .picture import Picture
from .dot_list import DotList
from .rectangle_list import RectangleList

ShapeType = _tp.Union[Dot, Rectangle, Picture]
