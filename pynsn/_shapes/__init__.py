import typing as _tp

from .dot import Dot
from .rectangle import Rectangle
from .picture import Picture
from .dot_list import DotList, BaseDotList
from .rectangle_list import RectangleList, RectangleLike, BaseRectangleList

ShapeType = _tp.Union[Dot, Rectangle, Picture]
