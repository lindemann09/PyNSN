from .. import arrays as _arrays
from . import fit as _adapt

from .fit import change_adapt_settings # make available


def average_diameter(dot_array, factor):
    if not isinstance(dot_array, _arrays.DotArray):
        raise TypeError("Scaling diameter is not possible for {}.".format(
            type(dot_array).__name__))
    if factor == 1:
        return dot_array
    value = dot_array.properties.average_dot_diameter * factor
    return _adapt.average_diameter(dot_array, value)

def average_rectangle_size(rect_array, factor):
    if not isinstance(rect_array, _arrays.RectangleArray):
        raise TypeError("Scaling rectangle size is not possible for {}.".format(
            type(rect_array).__name__))
    if factor == 1:
        return rect_array
    value = rect_array.properties.average_rectangle_size * factor
    return _adapt.average_rectangle_size(rect_array, value)


def total_surface_area(object_array, factor):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.total_surface_area * factor
    return _adapt.total_surface_area(object_array, value)


def field_area(object_array, factor, precision=None):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.field_area * factor
    return _adapt.field_area(object_array, value, precision=precision)

def coverage(object_array, factor,
             precision=None,
             FA2TA_ratio=None):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.converage * factor
    return _adapt.coverage(object_array, value,
                           precision=precision,
                           FA2TA_ratio=FA2TA_ratio)

def average_perimeter(object_array, factor):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.average_perimeter * factor
    return _adapt.average_perimeter(object_array, value)

def total_perimeter(object_array, factor):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.total_perimeter * factor
    return _adapt.total_perimeter(object_array, value)

def average_surface_area(object_array, factor):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.average_surface_area * factor
    return _adapt.average_surface_area(object_array, value)

def log_spacing(object_array, factor, precision=None):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.log_spacing * factor
    return _adapt.log_spacing(object_array, value, precision=precision)

def log_size(object_array, factor):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.log_size * factor
    return _adapt.log_size(object_array, value)

def sparcity(object_array, factor, precision=None):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.sparsity * factor
    return _adapt.sparcity(object_array, value, precision=precision)

def visual_property(object_array, feature, factor):
    if factor == 1:
        return object_array
    _arrays._check_object_array(object_array)
    value = object_array.properties.get(feature) * factor
    return _adapt.visual_property(object_array, property_flag=feature, value=value)
