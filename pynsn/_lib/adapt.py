import random as _random
import numpy as _np
from . import misc as _misc
from . import geometry as _geometry
from .arrays import DotArray as _DotArray
from .arrays import RectangleArray as _RectangleArray
from .arrays import _check_object_array
from .visual_features import VisualFeature as _VF
from . import adapt_settings as _settings

from .adapt_settings import change_adapt_settings # make available


#FIXME coverage for all
##FIXME FIELD_AREA_POSITIONS

def average_diameter(dot_array, value):
    # changes diameter
    if not isinstance(dot_array, _DotArray):
        raise TypeError("Adapting diameter is not possible for {}.".format(
            type(dot_array).__name__))
    scale = value / dot_array.features.average_dot_diameter
    dot_array._diameters = dot_array.diameters * scale
    dot_array.features.reset()
    return dot_array

def average_rectangle_size(rect_array, value):
    # changes diameter
    if not isinstance(rect_array, _RectangleArray):
        raise TypeError("Adapting rectangle size is not possible for {}.".format(
            type(rect_array).__name__))
    try:
        width, height = value
    except TypeError:
        raise TypeError("Value ({}) has to tuple of 2 numerical (width, height).".format(
            value))

    scale = _np.array([width / rect_array.features.average_rectangle_size[0],
                       height / rect_array.features.average_rectangle_size[1]])
    rect_array._sizes = rect_array._sizes * scale
    rect_array.features.reset()
    return rect_array

def total_surface_area(object_array, value):
    # changes diameter
    _check_object_array(object_array)
    a_scale = value / object_array.features.total_surface_area
    if isinstance(object_array, _DotArray):
        object_array._diameters = _np.sqrt(
            object_array.surface_areas * a_scale) * 2 / _np.sqrt(
                    _np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
    else: # rect
        object_array._sizes = object_array._sizes * _np.sqrt(a_scale)

    object_array.features.reset()
    return object_array

def field_area(object_array, value, precision=None):
    """changes the convex hull area to a desired size with certain precision

    uses scaling radial positions if field area has to be increased
    uses replacement of outer points (and later re-scaling)

    iterative method can takes some time.
    """

    _check_object_array(object_array)
    if precision is None:
        precision = _settings.DEFAULT_SPACING_PRECISION

    if object_array.features.field_area is _np.nan:
        return  object_array # not defined
    else:
        return _scale_field_area(object_array, value=value,
                                 precision=precision)

def _scale_field_area(object_array, value, precision):
    """change the convex hull area to a desired size by scale the polar
    positions  with certain precision

    iterative method can takes some time.

    Note: see doc string `_adapt_field_area`
    """
    _check_object_array(object_array)
    current = object_array.features.field_area

    if current is None:
        return  # not defined

    scale = 1  # find good scale
    step = 0.1
    if value < current:  # current too larger
        step *= -1

    # centered points
    old_center = object_array.center_of_mass()
    object_array._xy = object_array.xy - old_center
    centered_polar = _geometry.cartesian2polar(object_array.xy)

    # iteratively determine scale
    while abs(current - value) > precision:

        scale += step

        object_array._xy = _geometry.polar2cartesian(centered_polar * [scale, 1])
        object_array.features.reset()  # required at this point to recalc convex hull
        current = object_array.features.field_area

        if (current < value and step < 0) or \
                (current > value and step > 0):
            step *= -0.2  # change direction and finer grain

    object_array._xy = object_array.xy + old_center
    object_array.features.reset()
    return object_array

def coverage(object_array, value,
             precision=None,
             FA2TA_ratio=None):

    # FIXME check drifting outwards if extra space is small and adapt_FA2TA_ratio=1
    # FIXME when to realign, realignment changes field_area!
    """this function changes the area and remixes to get a desired density
    precision in percent between 1 < 0

    ratio_area_convex_hull_adaptation:
        ratio of adaptation via area or via convex_hull (between 0 and 1)

    """
    _check_object_array(object_array)

    print("WARNING: _adapt_coverage is a experimental ")
    # dens = convex_hull_area / total_surface_area
    if FA2TA_ratio is None:
        FA2TA_ratio = _settings.DEFAULT_ADAPT_FA2TA_RATIO
    elif FA2TA_ratio < 0 or FA2TA_ratio > 1:
        FA2TA_ratio = 0.5
    if precision is None:
        precision = _settings.DEFAULT_SPACING_PRECISION

    total_area_change100 = (value * object_array.features.field_area) - \
                           object_array.features.total_surface_area
    d_change_total_area = total_area_change100 * (1 - FA2TA_ratio)
    if abs(d_change_total_area) > 0:
        total_surface_area(object_array.features.total_surface_area + \
                           d_change_total_area)

    return field_area(object_array.features.total_surface_area / value,
               precision=precision)

def average_perimeter(object_array, value):
    _check_object_array(object_array)
    total_peri = value * object_array.features.numerosity
    return total_perimeter(object_array, total_peri)

def total_perimeter(object_array, value):
    if isinstance(object_array, _DotArray):
        tmp = value / (object_array.features.numerosity * _np.pi)
        return average_diameter(object_array, tmp)
    elif isinstance(object_array, _RectangleArray):
        scale = value / object_array.features.total_perimeter
        new_size = object_array.features.average_rectangle_size * scale
        return average_rectangle_size(object_array, new_size)
    else:
        _check_object_array(object_array)


def average_surface_area(object_array, value):
    _check_object_array(object_array)
    ta = object_array.features.numerosity * value
    return total_surface_area(object_array, ta)

def log_spacing(object_array, value, precision=None):
    _check_object_array(object_array)
    logfa = 0.5 * value + 0.5 * _misc.log2(
        object_array.features.numerosity)
    return field_area(object_array, value=2 ** logfa, precision=precision)

def log_size(object_array, value):
    _check_object_array(object_array)
    logtsa = 0.5 * value + 0.5 * _misc.log2(object_array.features.numerosity)
    return total_surface_area(object_array, 2 ** logtsa)

def sparcity(object_array, value, precision=None):
    _check_object_array(object_array)
    return field_area(object_array, value=value * object_array.features.numerosity,
                      precision=precision)

def visual_feature(object_array, feature, value):
    """
    adapt_properties: continuous property or list of continuous properties
    several properties to be adapted
    if adapt dot array is specified, array will be adapt to adapt_dot_array, otherwise
    the values defined in adapt_features is used.
    some adapting requires realignement to avoid overlaps. However,
    realigment might result in a different field area. Thus, realign after
    adapting for  Size parameter and realign before adapting space
    parameter.

    """

    # type check
    if not isinstance(feature, _VF):
        raise ValueError("{} is not a visual feature.".format(feature))

    # Adapt
    if feature == _VF.AV_DOT_DIAMETER:
        return average_diameter(object_array, value=value)

    elif feature == _VF.AV_RECT_SIZE:
        return average_rectangle_size(object_array, value=value)

    elif feature == _VF.AV_PERIMETER:
        return average_perimeter(object_array, value=value)

    elif feature == _VF.TOTAL_PERIMETER:
        return total_perimeter(object_array, value=value)

    elif feature == _VF.AV_SURFACE_AREA:
        return average_surface_area(object_array, value=value)

    elif feature == _VF.TOTAL_SURFACE_AREA:
        return total_surface_area(object_array, value=value)

    elif feature == _VF.LOG_SIZE:
        return log_size(object_array, value=value)

    elif feature == _VF.LOG_SPACING:
        return log_spacing(object_array, value=value)

    elif feature == _VF.SPARSITY:
        return sparcity(object_array, value=value)

    elif feature == _VF.FIELD_AREA:
        return field_area(object_array, value=value)

    elif feature == _VF.COVERAGE:
        return coverage(object_array, value=value)

    else:
        raise NotImplementedError("Not implemented for {}".format(
            feature.label()))