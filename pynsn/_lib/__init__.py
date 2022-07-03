from .base_array import BaseArray
from .dot_array import DotArray
from .rect_array import RectangleArray
from .shapes import Point, Dot, Rectangle
from .misc import PictureFile
from .misc import NoSolutionError


# helper for type checking and error raising error
def _check_base_array(obj):
    if not isinstance(obj, BaseArray):
        raise TypeError("DotArray, RectangleArray or AttributeArray expected, but not {}".format(
            type(obj).__name__))


def _check_object_array(obj):
    if not isinstance(obj, (DotArray, RectangleArray)):
        raise TypeError("DotArray or RectangleArray expected, but not {}".format(
            type(obj).__name__))

