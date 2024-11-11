import typing as _tp
from .shapes import Rectangle, Picture, PolygonShape, AbstractShape
from .circ_shapes import Point2D, Dot, Ellipse
from .colour import Colour


def dict_to_shape(d: _tp.Dict[str, _tp.Any]) -> _tp.Optional[AbstractShape]:
    """helper function converts dict to a shape.
    Returns None, if dict does not contain a shape"""

    if d["type"] == Dot.shape_type():
        return Dot.from_dict(d)
    elif d["type"] == Ellipse.shape_type():
        return Ellipse.from_dict(d)
    elif d["type"] == Rectangle.shape_type():
        return Rectangle.from_dict(d)
    elif d["type"] == PolygonShape.shape_type():
        return PolygonShape.from_dict(d)

    return None
