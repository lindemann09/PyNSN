from ._base_array import BaseArray
from ._dot_array import DotArray
from ._rect_array import RectangleArray

# helper for type checking and error raising error
def _check_base_array(obj):
    if not isinstance(obj, BaseArray):
        raise TypeError("DotArray, RectangleArray or AttributeArray expected, but not {}".format(
            type(obj).__name__))

def _check_object_array(obj):
    if not isinstance(obj, (DotArray, RectangleArray)):
        raise TypeError("DotArray or RectangleArray expected, but not {}".format(
            type(obj).__name__))

