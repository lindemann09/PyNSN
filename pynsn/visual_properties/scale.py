from . import fit as _fit
from .. import _lib

from .fit import change_fit_settings # make available


def average_diameter(dot_array, factor):
    if not isinstance(dot_array, _lib.DotArray):
        raise TypeError("Scaling diameter is not possible for {}.".format(
            type(dot_array).__name__))
    if factor == 1:
        return dot_array
    value = dot_array.properties.average_dot_diameter * factor
    return _fit.average_diameter(dot_array, value)

def average_rectangle_size(rect_array, factor):
    if not isinstance(rect_array, _lib.RectangleArray):
        raise TypeError("Scaling rectangle size is not possible for {}.".format(
            type(rect_array).__name__))
    if factor == 1:
        return rect_array
    value = rect_array.properties.average_rectangle_size * factor
    return _fit.average_rectangle_size(rect_array, value)


def total_surface_area(object_array, factor):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.total_surface_area * factor
    return _fit.total_surface_area(object_array, value)


def field_area(object_array, factor, precision=None):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.field_area * factor
    return _fit.field_area(object_array, value, precision=precision)

def coverage(object_array, factor,
             precision=None,
             FA2TA_ratio=None):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.converage * factor
    return _fit.coverage(object_array, value,
                        precision=precision,
                        FA2TA_ratio=FA2TA_ratio)

def average_perimeter(object_array, factor):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.average_perimeter * factor
    return _fit.average_perimeter(object_array, value)

def total_perimeter(object_array, factor):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.total_perimeter * factor
    return _fit.total_perimeter(object_array, value)

def average_surface_area(object_array, factor):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.average_surface_area * factor
    return _fit.average_surface_area(object_array, value)

def log_spacing(object_array, factor, precision=None):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.log_spacing * factor
    return _fit.log_spacing(object_array, value, precision=precision)

def log_size(object_array, factor):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.log_size * factor
    return _fit.log_size(object_array, value)

def sparcity(object_array, factor, precision=None):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.sparsity * factor
    return _fit.sparcity(object_array, value, precision=precision)

def visual_property(object_array, feature, factor):
    if factor == 1:
        return object_array
    _lib._check_object_array(object_array)
    value = object_array.properties.get(feature) * factor
    return _fit.visual_property(object_array, property_flag=feature, value=value)
