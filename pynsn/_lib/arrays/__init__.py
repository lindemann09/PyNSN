from ._generic_object_array import GenericObjectArray
from ._dot_array import DotArray
from ._rect_array import RectangleArray

# helper functions
def _check_generic_array(obj):
    if not isinstance(obj, GenericObjectArray):
        raise TypeError("DotArray, RectangleArray or GenericObjectArray expected, but not {}".format(
                        type(obj).__name__))

def _check_object_array(obj):
    if not isinstance(obj, (DotArray, RectangleArray)):
        raise TypeError("DotArray or RectangleArray expected, but not {}".format(
                        type(obj).__name__))
