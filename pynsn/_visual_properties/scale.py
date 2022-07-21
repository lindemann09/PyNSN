from __future__ import annotations

from . import fit as _fit
from .. import _arrays
from .._lib import lib_typing as _tp

from .properties import VisualPropertyFlag as flags

def average_diameter(dot_array: _fit._ObjectArrayType, factor: float) -> None:
    if not isinstance(dot_array, _arrays.DotArray):
        raise TypeError("Scaling diameter is not possible for {}.".format(
            type(dot_array).__name__))
    if factor == 1:
        return
    value = dot_array.properties.average_dot_diameter * factor
    _fit.average_diameter(dot_array, value)


def average_rectangle_size(rect_array: _fit._ObjectArrayType, factor: float) -> None:
    if not isinstance(rect_array, _arrays.RectangleArray):
        raise TypeError("Scaling rectangle size is not possible for {}.".format(
            type(rect_array).__name__))
    if factor == 1:
        return
    value = rect_array.properties.average_rectangle_size * factor
    _fit.average_rectangle_size(rect_array, value)


def total_surface_area(object_array: _fit._ObjectArrayType, factor: float) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.total_surface_area * factor
    _fit.total_surface_area(object_array, value)


def field_area(object_array: _fit._ObjectArrayType, factor: float,
               precision: _tp.OptFloat = None) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.field_area * factor
    _fit.field_area(object_array, value, precision=precision)


def coverage(object_array: _fit._ObjectArrayType, factor: float,
             precision: _tp.OptFloat = None,
             FA2TA_ratio: _tp.OptFloat = None) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.converage * factor
    _fit.coverage(object_array, value,
                  precision=precision,
                  FA2TA_ratio=FA2TA_ratio)


def average_perimeter(object_array: _fit._ObjectArrayType, factor: float) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.average_perimeter * factor
    _fit.average_perimeter(object_array, value)


def total_perimeter(object_array: _fit._ObjectArrayType, factor: float) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.total_perimeter * factor
    _fit.total_perimeter(object_array, value)


def average_surface_area(object_array: _fit._ObjectArrayType, factor: float) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.average_surface_area * factor
    _fit.average_surface_area(object_array, value)


def log_spacing(object_array: _fit._ObjectArrayType, factor: float,
                precision: _tp.OptFloat = None) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.log_spacing * factor
    _fit.log_spacing(object_array, value, precision=precision)


def log_size(object_array: _fit._ObjectArrayType, factor: float) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.log_size * factor
    _fit.log_size(object_array, value)


def sparcity(object_array: _fit._ObjectArrayType, factor: float,
             precision: _tp.OptFloat = None) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.sparsity * factor
    _fit.sparcity(object_array, value, precision=precision)


def visual_property(object_array: _fit._ObjectArrayType, feature: flags, factor: float) -> None:
    if factor == 1:
        return
    _fit._check_object_array(object_array)
    value = object_array.properties.get(feature) * factor
    _fit.visual_property(object_array, property_flag=feature, value=value)
