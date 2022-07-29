import typing as _tp

from .point import Point
from .dot import Dot
from .rectangle import Rectangle
from .picture_file import PictureFile
ShapeType = _tp.Union[Point, Rectangle, Dot]
