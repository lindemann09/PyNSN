# FIXME not yet used implemented for Rectangles
import random as _random
import numpy as _np
from ._lib import misc as _misc
from ._lib import geometry as _geometry
from ._lib.arrays import DotArray as _DotArray
from ._lib.arrays import RectangleArray as _RectangleArray
from ._lib.visual_features import VisualFeatures as _VF

ITERATIVE_CONVEX_HULL_MODIFICATION = False
TAKE_RANDOM_DOT_FROM_CONVEXHULL = True  # TODO maybe method for modification, set via GUI
DEFAULT_SPACING_PRECISION = 0.0001
DEFAULT_MATCH_FA2TA_RATIO = 0.5

def change_settings(iterative_convex_hull_modification=None,
                    take_random_dot_from_convexhull=None,
                    default_spacing_precision=None,
                    default_match_fa2ta_ratio=None):
    """Changing class settings of feature matcher.

    This changes the settings of all feature matcher.


    Parameters
    ----------
    iterative_convex_hull_modification
    take_random_dot_from_convexhull
    default_spacing_precision
    default_match_fa2ta_ratio

    Returns
    -------

    """
    global ITERATIVE_CONVEX_HULL_MODIFICATION
    global TAKE_RANDOM_DOT_FROM_CONVEXHULL
    global DEFAULT_MATCH_FA2TA_RATIO
    global DEFAULT_SPACING_PRECISION
    if isinstance(iterative_convex_hull_modification, bool):
        ITERATIVE_CONVEX_HULL_MODIFICATION = iterative_convex_hull_modification
    if isinstance(take_random_dot_from_convexhull, bool):
        TAKE_RANDOM_DOT_FROM_CONVEXHULL = take_random_dot_from_convexhull
    if isinstance(default_spacing_precision, float):
        DEFAULT_SPACING_PRECISION = default_spacing_precision
    if isinstance(default_match_fa2ta_ratio, float):
        DEFAULT_MATCH_FA2TA_RATIO = default_match_fa2ta_ratio


def item_diameter(dot_array, value):
    # changes diameter
    assert isinstance(dot_array, _DotArray)
    scale = value / dot_array.features.mean_item_diameter
    dot_array._diameters = dot_array.diameters * scale
    dot_array.features.reset()
    return dot_array

def total_surface_area(object_array, value):
    # changes diameter
    assert isinstance(object_array, (_DotArray, _RectangleArray))

    a_scale = (value / object_array.features.total_surface_area)
    if isinstance(object_array, _DotArray):
        object_array._diameters = _np.sqrt(
            object_array.surface_areas * a_scale) * 2 / _np.sqrt(
                    _np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
    object_array.features.reset()
    return object_array

def field_area(object_array, value, precision=None, use_scaling_only=False):
    """changes the convex hull area to a desired size with certain precision

    uses scaling radial positions if field area has to be increased
    uses replacement of outer points (and later re-scaling)

    iterative method can takes some time.
    """

    # PROCEDURE
    #
    # increasing field area:
    #   * merely scales polar coordinates of Dots
    #   * uses __scale_field_area
    #
    # decreasing field area:
    #   a) iterative convex hull modification
    #      1. iteratively replacing outer dots to the side (random  pos.)
    #         (resulting FA is likely to small)
    #      2. increase FA by scaling to match precisely
    #         inside the field area
    #      - this methods results in very angular dot arrays, because it
    #           prefers a solution with a small number of convex hull
    #           dots
    #   b) eccentricity criterion
    #      1. determining circle with the required field area
    #      2. replacing all dots outside this circle to the inside
    #         (random pos.) (resulting FA is likely to small)
    #      3. increase FA by scaling to match precisely
    #      - this method will result is rather circular areas

    assert isinstance(object_array, _DotArray)
    if precision is None:
        precision = DEFAULT_SPACING_PRECISION

    if object_array.features.field_area is None:
        return  object_array # not defined
    elif value > object_array.features.field_area or use_scaling_only:
        # field area is currently too small or scaling is enforced
        return _scale_field_area(object_array, value=value,
                                 precision=precision)
    elif value < object_array.features.field_area:
        # field area is too large
        _decrease_field_area_by_replacement(object_array,
                                            max_field_area=value,
                                            iterative_convex_hull_modification=
            ITERATIVE_CONVEX_HULL_MODIFICATION)
        # ..and rescaling to avoid to compensate for possible too
        # strong decrease
        return _scale_field_area(object_array, value=value, precision=precision)
    else:
        return object_array

def _decrease_field_area_by_replacement(object_array, max_field_area,
                                        iterative_convex_hull_modification):
    """decreases filed area by recursively moving the most outer point
    to some more central free position (avoids overlapping)

    return False if not possible else True

    Note: see doc string `_match_field_area`

    """
    assert isinstance(object_array, _DotArray)
    # centered points
    old_center = object_array.center_of_outer_positions
    object_array._xy = object_array._xy - old_center

    removed_dots = []

    if iterative_convex_hull_modification:
        while object_array.features.field_area > max_field_area:
            # remove one random outer dot and remember it
            indices = object_array.features.convex_hull.indices
            if not TAKE_RANDOM_DOT_FROM_CONVEXHULL:
                # most outer dot from convex hull
                radii_outer_dots = _geometry.cartesian2polar(
                                object_array.xy[indices],
                                radii_only=True)
                i = _np.where(radii_outer_dots == max(radii_outer_dots))[0]
                idx = indices[i][0]
            else:
                # remove random
                idx = indices[_random.randint(0, len(indices) - 1)]

            removed_dots.extend(object_array.get(indices=[idx]))
            object_array.delete(idx)

        # add dots to free pos inside the convex hall
        for d in removed_dots:
            try:
                d._xy = object_array.random_free_position(d.diameter,
                                                          allow_overlapping=False,
                                                          prefer_inside_field_area=True)
            except StopIteration as e:
                raise StopIteration("Can't find a free position while decreasing field area.\n" +
                                   "n={}; current FA={}, max_FA={}".format(
                                       object_array.features.numerosity + 1,
                                       object_array.features.field_area,
                                       max_field_area))

            object_array.add(d)

    else:
        # eccentricity criterion
        max_radius = _np.sqrt(max_field_area / _np.pi)  # for circle with
        # required FA
        idx = _np.where(_geometry.cartesian2polar(object_array.xy, radii_only=True) > max_radius)[0]
        removed_dots.extend(object_array.get(indices=idx))
        object_array.delete(idx)

        # add inside the circle
        min_dist = object_array.target_array_radius - max_radius + 1
        for d in removed_dots:
            try:
                d._xy = object_array.random_free_position(d.diameter,
                                                          prefer_inside_field_area=False,
                                                          allow_overlapping=False,
                                                          min_distance_area_boarder=min_dist)
            except StopIteration as e:
                raise StopIteration(
                    "Can't find a free position while decreasing field area.\n" + \
                    "n={}; current FA={}, max_FA={}".format(
                        object_array.features.numerosity + 1,
                        object_array.features.field_area,
                        max_field_area))
            object_array.add(d)

    object_array._xy = object_array.xy + old_center
    object_array.features.reset()
    return object_array

def _scale_field_area(object_array, value, precision):
    """change the convex hull area to a desired size by scale the polar
    positions  with certain precision

    iterative method can takes some time.

    Note: see doc string `_match_field_area`
    """
    assert isinstance(object_array, _DotArray)
    current = object_array.features.field_area

    if current is None:
        return  # not defined

    scale = 1  # find good scale
    step = 0.1
    if value < current:  # current too larger
        step *= -1

    # centered points
    old_center = object_array.center_of_outer_positions
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
             match_FA2TA_ratio=None):

    # FIXME check drifting outwards if extra space is small and match_FA2TA_ratio=1
    # FIXME when to realign, realignment changes field_area!
    """this function changes the area and remixes to get a desired density
    precision in percent between 1 < 0

    ratio_area_convex_hull_adaptation:
        ratio of adaptation via area or via convex_hull (between 0 and 1)

    """
    assert isinstance(object_array, _DotArray)

    print("WARNING: _match_coverage is a experimental ")
    # dens = convex_hull_area / total_surface_area
    if match_FA2TA_ratio is None:
        match_FA2TA_ratio = DEFAULT_MATCH_FA2TA_RATIO
    elif match_FA2TA_ratio < 0 or match_FA2TA_ratio > 1:
        match_FA2TA_ratio = 0.5
    if precision is None:
        precision = DEFAULT_SPACING_PRECISION

    total_area_change100 = (value * object_array.features.field_area) - \
                           object_array.features.total_surface_area
    d_change_total_area = total_area_change100 * (1 - match_FA2TA_ratio)
    if abs(d_change_total_area) > 0:
        total_surface_area(object_array.features.total_surface_area + \
                           d_change_total_area)

    field_area(object_array.features.total_surface_area / value,
               precision=precision)

def item_perimeter(object_array, value):
    return item_diameter(object_array, value / _np.pi)

def total_perimeter(object_array, value):
    assert isinstance(object_array, _DotArray)
    tmp = value / (object_array.features.numerosity * _np.pi)
    return item_diameter(object_array, tmp)

def item_surface_area(object_array, value):
    assert isinstance(object_array, _DotArray)
    ta = object_array.features.numerosity * value
    return total_surface_area(object_array, ta)

def log_spacing(object_array, value, precision=None):
    assert isinstance(object_array, _DotArray)
    logfa = 0.5 * value + 0.5 * _misc.log2(
        object_array.features.numerosity)
    return field_area(object_array, value=2 ** logfa, precision=precision)

def log_size(object_array, value):
    assert isinstance(object_array, _DotArray)
    logtsa = 0.5 * value + 0.5 * _misc.log2(object_array.features.numerosity)
    return total_surface_area(object_array, 2 ** logtsa)

def sparcity(object_array, value, precision=None):
    assert isinstance(object_array, _DotArray)
    return field_area(object_array, value=value * object_array.features.numerosity,
                      precision=precision)

def visual_feature(object_array, feature, value):
    """
    match_properties: continuous property or list of continuous properties
    several properties to be matched
    if match dot array is specified, array will be match to match_dot_array, otherwise
    the values defined in match_features is used.
    some matching requires realignement to avoid overlaps. However,
    realigment might result in a different field area. Thus, realign after
    matching for  Size parameter and realign before matching space
    parameter.

    """

    # type check
    if feature not in _VF.ALL_FEATURES:
        raise ValueError("{} is not a visual feature.".format(feature))

    # Adapt
    if feature == _VF.ITEM_DIAMETER:
        return item_diameter(object_array, value=value)

    elif feature == _VF.ITEM_PERIMETER:
        return item_perimeter(object_array, value=value)

    elif feature == _VF.TOTAL_PERIMETER:
        return total_perimeter(object_array, value=value)

    elif feature == _VF.ITEM_SURFACE_AREA:
        return item_surface_area(object_array, value=value)

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
