from .arrays import DotArray as _DotArray
from .arrays import RectangleArray as _RectangleArray
from . import adapt as _adapt

from .adapt_settings import change_adapt_settings # make available


def average_diameter(dot_array, factor):
    if not isinstance(dot_array, _DotArray):
        raise TypeError("Scaling diameter is not possible for {}.".format(
            type(dot_array).__name__))
    value = dot_array.features.average_dot_diameter * factor
    return _adapt.average_diameter(dot_array, value)

def average_rect_size(rect_array, factor):
    if not isinstance(rect_array, _RectangleArray): # FIXME make check type functions
        raise TypeError("Scaling rectangle size is not possible for {}.".format(
            type(rect_array).__name__))
    value = rect_array.features.average_rectangle_size * factor
    return _adapt.average_rect_size(rect_array, value)


def total_surface_area(object_array, factor):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.total_surface_area * factor
    return _adapt.total_surface_area(object_array, value)


def field_area(object_array, factor, precision=None, use_scaling_only=False):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.field_area * factor
    return _adapt.field_area(object_array, value)

def coverage(object_array, factor,
             precision=None,
             FA2TA_ratio=None):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.converage * factor
    return _adapt.coverage(object_array, value,
                           precision=precision,
                           FA2TA_ratio=FA2TA_ratio)

def average_perimeter(object_array, factor):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.average_perimeter * factor
    return _adapt.average_perimeter(object_array, value)

def total_perimeter(object_array, factor):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.total_perimeter * factor
    return _adapt.total_perimeter(object_array, value)

def average_surface_area(object_array, factor):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.average_surface_area * factor
    return _adapt.average_surface_area(object_array, value)

def log_spacing(object_array, factor, precision=None):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.log_spacing * factor
    return _adapt.log_spacing(object_array, value, precision=precision)

def log_size(object_array, factor):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.log_size * factor
    return _adapt.log_size(object_array, value)

def sparcity(object_array, factor, precision=None):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.sparsity * factor
    return _adapt.sparcity(object_array, value, precision=precision)

def visual_feature(object_array, feature, factor):
    assert isinstance(object_array, (_DotArray, _RectangleArray))
    value = object_array.features.get(feature) * factor
    return _adapt.visual_feature(object_array, feature=feature, value=value)
