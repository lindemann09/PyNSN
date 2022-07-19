import numpy as _np
from math import log2 as _log2
from .._lib import geometry as _geometry
from .._lib import rng as _rng
from .. import _lib
from .. import constants
from ..exceptions import NoSolutionError as _NoSolutionError
from ._properties import VisualPropertyFlag as _flags


def change_fit_settings(default_spacing_precision=None,
                        default_fa2ta_ratio=None):
    """Changing class settings of property fitting.

    This changes the settings of the property fitting.


    Parameters
    ----------
    default_spacing_precision
    default_fa2ta_ratio

    Returns
    -------

    """

    if default_spacing_precision is not None:
        constants.DEFAULT_FIT_SPACING_PRECISION = float(default_spacing_precision)
    if default_fa2ta_ratio is not None:
        constants.DEFAULT_FIT_FA2TA_RATIO = float(default_fa2ta_ratio)

#FIXME coverage for all


def numerosity(object_array, value, keep_field_area=False):
    """

    """
    _lib._check_object_array(object_array)

    # make a copy for the deviant
    if value <= 0:
        object_array.clear()
    else:
        # add or remove random dots
        change_numerosity = value - object_array._properties.numerosity
        for _ in range(abs(change_numerosity)):
            if change_numerosity < 0:
                # remove dots
                if keep_field_area:
                    # find a random object that is not in convex hull
                    delete_id = None
                    ch = object_array._properties.convex_hull.indices
                    rnd_seq = list(range(0, object_array._properties.numerosity))
                    _rng.generator.shuffle(rnd_seq)
                    for x in rnd_seq:
                        if x not in ch:
                            delete_id = x
                            break
                    if delete_id is None:
                        raise _NoSolutionError("Can't increase numeroisty, while keeping field area.")
                else:
                    delete_id = _rng.generator.integers(0, object_array._properties.numerosity)

                object_array.delete(delete_id)

            else:
                # add dot: copy a random dot
                clone_id = _rng.generator.integers(0, object_array._properties.numerosity)
                rnd_object = next(object_array.iter_objects(clone_id))
                try:
                    rnd_object = object_array.get_free_position(
                        ref_object=rnd_object, allow_overlapping=False,
                        inside_convex_hull=keep_field_area)
                except _NoSolutionError:
                    # no free position
                    raise _NoSolutionError("Can't increase numerosity. No free position found.")

                object_array.add([rnd_object])

    return object_array


def average_diameter(dot_array, value):
    # changes diameter
    if not isinstance(dot_array, _lib.DotArray):
        raise TypeError("Adapting diameter is not possible for {}.".format(
            type(dot_array).__name__))
    scale = value / dot_array.properties.average_dot_diameter
    dot_array._diameter = dot_array.diameter * scale
    dot_array.properties.reset()
    return dot_array


def average_rectangle_size(rect_array, value):
    # changes diameter
    if not isinstance(rect_array, _lib.RectangleArray):
        raise TypeError("Adapting rectangle size is not possible for {}.".format(
            type(rect_array).__name__))
    try:
        width, height = value
    except TypeError:
        raise TypeError("Value ({}) has to tuple of 2 numerical (width, height).".format(
            value))

    scale = _np.array([width / rect_array.properties.average_rectangle_size[0],
                       height / rect_array.properties.average_rectangle_size[1]])
    rect_array._sizes = rect_array._sizes * scale
    rect_array.properties.reset()
    return rect_array


def total_surface_area(object_array, value):
    # changes diameter
    _lib._check_object_array(object_array)
    a_scale = value / object_array.properties.total_surface_area
    if isinstance(object_array, _lib.DotArray):
        object_array._diameter = _np.sqrt(
            object_array.surface_areas * a_scale) * 2 / _np.sqrt(
                    _np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
    else: # rect
        object_array._sizes = object_array._sizes * _np.sqrt(a_scale)

    object_array.properties.reset()
    return object_array


def field_area(object_array, value, precision=None):
    """changes the convex hull area to a desired size with certain precision

    uses scaling radial positions if field area has to be increased
    uses replacement of outer points (and later re-scaling)

    iterative method can takes some time.
    """

    _lib._check_object_array(object_array)
    if precision is None:
        precision = constants.DEFAULT_FIT_SPACING_PRECISION

    if object_array.properties.field_area is _np.nan:
        return  object_array # not defined
    else:
        return _scale_field_area(object_array, value=value,
                                 precision=precision)


def _scale_field_area(object_array, value, precision):
    """change the convex hull area to a desired size by scale the polar
    positions with certain precision

    iterative method can takes some time.

    Note: see doc string `field_area`
    """
    _lib._check_object_array(object_array)
    current = object_array.properties.field_area

    if current is None:
        return  # not defined

    scale = 1  # find good scale
    step = 0.1
    if value < current:  # current too larger
        step *= -1

    # centered points
    old_center = object_array.get_
    object_array._xy = object_array.xy - old_center
    centered_polar = _geometry.cartesian2polar(object_array.xy)

    # iteratively determine scale
    while abs(current - value) > precision:

        scale += step

        object_array._xy = _geometry.polar2cartesian(centered_polar * [scale, 1])
        object_array.properties.reset()  # required at this point to recalc convex hull
        current = object_array.properties.field_area

        if (current < value and step < 0) or \
                (current > value and step > 0):
            step *= -0.2  # change direction and finer grain

    object_array._xy = object_array.xy + old_center
    object_array.properties.reset()
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
    _lib._check_object_array(object_array)

    print("WARNING: _adapt_coverage is a experimental ")
    # dens = convex_hull_area / total_surface_area
    if FA2TA_ratio is None:
        FA2TA_ratio = constants.DEFAULT_FIT_FA2TA_RATIO
    elif FA2TA_ratio < 0 or FA2TA_ratio > 1:
        FA2TA_ratio = 0.5
    if precision is None:
        precision = constants.DEFAULT_FIT_SPACING_PRECISION

    total_area_change100 = (value * object_array.properties.field_area) - \
                           object_array.properties.total_surface_area
    d_change_total_area = total_area_change100 * (1 - FA2TA_ratio)
    if abs(d_change_total_area) > 0:
        total_surface_area(object_array.properties.total_surface_area + \
                           d_change_total_area)

    return field_area(object_array.properties.total_surface_area / value,
               precision=precision)


def average_perimeter(object_array, value):
    _lib._check_object_array(object_array)
    total_peri = value * object_array.properties.numerosity
    return total_perimeter(object_array, total_peri)


def total_perimeter(object_array, value):
    if isinstance(object_array, _lib.DotArray):
        tmp = value / (object_array.properties.numerosity * _np.pi)
        return average_diameter(object_array, tmp)
    elif isinstance(object_array, _lib.RectangleArray):
        scale = value / object_array.properties.total_perimeter
        new_size = object_array.properties.average_rectangle_size * scale
        return average_rectangle_size(object_array, new_size)
    else:
        _lib._check_object_array(object_array)


def average_surface_area(object_array, value):
    _lib._check_object_array(object_array)
    ta = object_array.properties.numerosity * value
    return total_surface_area(object_array, ta)


def log_spacing(object_array, value, precision=None):
    _lib._check_object_array(object_array)
    logfa = 0.5 * value + 0.5 * _log2(
        object_array.properties.numerosity)
    return field_area(object_array, value=2 ** logfa, precision=precision)


def log_size(object_array, value):
    _lib._check_object_array(object_array)
    logtsa = 0.5 * value + 0.5 * _log2(object_array.properties.numerosity)
    return total_surface_area(object_array, 2 ** logtsa)


def sparcity(object_array, value, precision=None):
    _lib._check_object_array(object_array)
    return field_area(object_array, value=value * object_array.properties.numerosity,
                      precision=precision)


def visual_property(object_array, property_flag, value):
    """
    adapt_properties: continuous property or list of continuous properties
    several properties to be adapted
    if adapt dot array is specified, array will be adapt to adapt_dot_array, otherwise
    the values defined in adapt_properties is used.
    some adapting requires realignement to avoid overlaps. However,
    realigment might result in a different field area. Thus, realign after
    adapting for  Size parameter and realign before adapting space
    parameter.

    """

    # type check
    if not isinstance(property_flag, _flags):
        raise ValueError("{} is not a visual feature.".format(property_flag))

    # Adapt
    if property_flag == _flags.AV_DOT_DIAMETER:
        return average_diameter(object_array, value=value)

    elif property_flag == _flags.NUMEROSITY:
        return numerosity(object_array, value=value)

    elif property_flag == _flags.AV_PERIMETER:
        return average_perimeter(object_array, value=value)

    elif property_flag == _flags.TOTAL_PERIMETER:
        return total_perimeter(object_array, value=value)

    elif property_flag == _flags.AV_SURFACE_AREA:
        return average_surface_area(object_array, value=value)

    elif property_flag == _flags.TOTAL_SURFACE_AREA:
        return total_surface_area(object_array, value=value)

    elif property_flag == _flags.LOG_SIZE:
        return log_size(object_array, value=value)

    elif property_flag == _flags.LOG_SPACING:
        return log_spacing(object_array, value=value)

    elif property_flag == _flags.SPARSITY:
        return sparcity(object_array, value=value)

    elif property_flag == _flags.FIELD_AREA:
        return field_area(object_array, value=value)

    elif property_flag == _flags.COVERAGE:
        return coverage(object_array, value=value)

    else:
        raise NotImplementedError("Not implemented for {}".format(
            property_flag.label()))

# TODO "visual test" (eye inspection) of fitting rect arrays
