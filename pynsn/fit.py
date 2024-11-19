"""Fitting module"""
import typing as _tp

import numpy as _np
import shapely as _shp

from . import NSNStimulus as _NSNStimulus
from . import NSNStimulusPair as _NSNStimulusPair
from . import VisProp as _VisProp
from . import _misc
from . import defaults as _defaults
from .exceptions import NoSolutionError as _NoSolutionError
from .rnd import generator as _rnd_generator


def total_surface_area(stim: _NSNStimulus,
                       value: _tp.Union[float, _np.float64]) -> None:
    """Set surface area.

    Resize all object to fit a specific surface area

    Args:
        value: surface area
    """
    stim.sizes = stim.sizes * \
        _np.sqrt((value / stim.properties.total_surface_area))


def total_perimeter(stim: _NSNStimulus,
                    value: _tp.Union[float, _np.float64]) -> None:
    """fit the total parameter of the stimulus"""
    stim.sizes = stim.sizes * (value / stim.properties.total_perimeter)


def average_perimeter(stim: _NSNStimulus, value: _tp.Union[float, _np.float64]) -> None:
    """fit the average parameter of the stimulus"""
    total_perimeter(stim, value * stim.n_shapes)


def average_surface_area(stim: _NSNStimulus, value: _tp.Union[float, _np.float64]) -> None:
    """fits the average surface area of the stimulus"""
    total_surface_area(stim, stim.n_shapes * value)


def numerosity(stim: _NSNStimulus, value: int,
               keep_convex_hull: bool = False,
               max_iterations: _tp.Optional[int] = None) -> None:
    """
    fitting the numerosity
    """

    # make a copy for the deviant
    if value <= 0:
        stim.shapes_clear()
    else:
        # add or remove random dots
        change_numerosity = value - stim.n_shapes
        if keep_convex_hull and change_numerosity < 0:
            # find shapes touching the convex hull (ch_shapes)
            ring = _shp.get_exterior_ring(stim.convex_hull.polygon)
            ch_shapes = _np.flatnonzero(_shp.intersects(
                stim.polygons, ring))
        else:
            ch_shapes = None

        for _ in range(abs(change_numerosity)):
            if change_numerosity < 0:
                # remove dots
                if ch_shapes is not None:
                    # find a random object that is not in convex hull
                    rnd_seq = _np.arange(stim.n_shapes)
                    _rnd_generator.shuffle(rnd_seq)
                    delete_id = None
                    for x in rnd_seq:
                        if x not in ch_shapes:
                            delete_id = x
                            break
                    if delete_id is None:
                        raise _NoSolutionError(
                            "Can't increase numerosity, while keeping field area.")
                else:
                    delete_id = _rnd_generator.integers(0, stim.n_shapes)

                stim.shape_delete(delete_id)

            else:
                # add dot: copy a random dot
                clone_id = _rnd_generator.integers(0, stim.n_shapes)
                rnd_object = stim.shapes[clone_id]
                try:
                    rnd_object = stim.random_free_position(
                        shape=rnd_object, ignore_overlaps=False,
                        inside_convex_hull=keep_convex_hull,
                        max_iterations=max_iterations)
                except _NoSolutionError as err:
                    # no free position
                    raise _NoSolutionError(
                        "Can't increase numerosity. No free position found.") from err

                stim.shape_add(rnd_object)


def field_area(stim: _NSNStimulus,
               value: _tp.Union[float, _np.float64],
               precision: _tp.Optional[_tp.Union[float, _np.float64]] = None) -> None:
    """changes the convex hull area to a desired size with certain precision

    uses scaling radial positions if field area has to be increased
    uses replacement of outer points (and later re-scaling)

    iterative method can takes some time.
    """

    if precision is None or _np.isnan(precision):
        precision = _defaults.FIT_SPACING_PRECISION

    current = stim.convex_hull.area
    if stim.n_shapes < 3 or current == 0:
        return None  # not defined

    scale = 1  # find good scale
    step = 0.1
    if value < current:  # current too larger
        step *= -1

    # centred points
    old_center = stim.convex_hull.centroid
    stim.xy = stim.xy - old_center
    centered_polar = _misc.cartesian2polar(stim.xy)

    # iteratively determine scale
    while abs(current - value) > precision:
        scale += step

        stim.xy = _misc.polar2cartesian(centered_polar * [scale, 1])
        current = stim.convex_hull.area

        if (current < value and step < 0) or (current > value and step > 0):
            step *= -0.2  # change direction and finer grain

    stim.xy = stim.xy + old_center  # move back


def coverage(stim: _NSNStimulus, value: _tp.Union[float, _np.float64],
             precision: _tp.Optional[float] = None,
             fa2ta_ratio: _tp.Optional[float] = None) -> None:
    """

    Parameters
    ----------
    value
    precision
    fa2ta_ratio

    Returns
    -------

    """

    # TODO check drifting outwards if extra space is small and adapt_FA2TA_ratio=1
    # when to realign, realignment changes field_area!
    # """this function changes the area and remixes to get a desired density
    # precision in percent between 1 < 0
    #
    # ratio_area_convex_hull_adaptation:
    #    ratio of adaptation via area or via convex_hull (between 0 and 1)

    print("WARNING: _adapt_coverage is a experimental ")
    # dens = convex_hull_area / total_surface_area
    if fa2ta_ratio is None:
        fa2ta_ratio = _defaults.FIT_FA2TA_RATIO
    elif fa2ta_ratio < 0 or fa2ta_ratio > 1:
        fa2ta_ratio = 0.5
    if precision is None:
        precision = _defaults.FIT_SPACING_PRECISION

    ta = stim.properties.total_surface_area  # total area
    ta_change100 = (value * stim.properties.field_area) - ta
    d_ta_change = ta_change100 * (1 - fa2ta_ratio)
    if abs(d_ta_change) > 0:
        total_surface_area(stim, ta + d_ta_change)

    field_area(stim,
               value=stim.properties.total_surface_area / value,
               precision=precision)


def log_spacing(stim: _NSNStimulus,
                value: _tp.Union[float, _np.float64],
                precision: _tp.Optional[float] = None) -> None:
    """

    Parameters
    ----------
    value
    precision

    Returns
    -------

    """
    log_fa = 0.5 * value + 0.5 * _np.log2(stim.n_shapes)
    field_area(stim, value=2**log_fa, precision=precision)


def log_size(stim: _NSNStimulus,
             value: _tp.Union[float, _np.float64]) -> None:
    """

    Parameters
    ----------
    value

    Returns
    -------

    """
    log_tsa = 0.5 * value + 0.5 * _np.log2(stim.n_shapes)
    total_surface_area(stim, value=2**log_tsa)


def sparsity(stim: _NSNStimulus,
             value: _tp.Union[float, _np.float64], precision=None) -> None:
    """

    Parameters
    ----------
    value
    precision

    Returns
    -------

    """
    field_area(stim, value=value * stim.n_shapes, precision=precision)


def property_adapt(stim: _NSNStimulus,
                   prop: _VisProp,
                   value: _tp.Union[float, _np.float64]) -> _tp.Any:
    """Adapt visual property `prop` of NSNStimulus


    Note:
    Realignment might be required after adapting visual properties
    """

    if not isinstance(prop, _VisProp):
        raise ValueError(
            f"{prop} is not a visual feature.")

    # Adapt
    if prop == _VisProp.NUMEROSITY:
        return numerosity(stim, value=int(value))

    elif prop == _VisProp.AV_PERIMETER:
        return average_perimeter(stim, value=value)

    elif prop == _VisProp.TOTAL_PERIMETER:
        return total_perimeter(stim, value=value)

    elif prop == _VisProp.AV_SURFACE_AREA:
        return average_surface_area(stim, value=value)

    elif prop == _VisProp.TOTAL_SURFACE_AREA:
        return total_surface_area(stim, value=value)

    elif prop == _VisProp.LOG_SIZE:
        return log_size(stim, value=value)

    elif prop == _VisProp.LOG_SPACING:
        return log_spacing(stim, value=value)

    elif prop == _VisProp.SPARSITY:
        return sparsity(stim, value=value)

    elif prop == _VisProp.FIELD_AREA:
        return field_area(stim, value=value)

    elif prop == _VisProp.COVERAGE:
        return coverage(stim, value=value)
    else:
        raise NotImplementedError(
            f"Not implemented for {str(prop)}"
        )


def property_match(stim: _NSNStimulus,
                   ref: _NSNStimulus,
                   prop: _VisProp) -> _tp.Any:
    """
    Match the visual property of stimulus `stim` to stimulus `ref`
    """
    # type check
    if not isinstance(prop, _VisProp):
        raise ValueError(
            f"{prop} is not a visual feature.")

    rp = ref.properties
    # Adapt
    if prop == _VisProp.NUMEROSITY:
        return numerosity(stim, value=rp.numerosity)

    elif prop == _VisProp.AV_PERIMETER:
        return average_perimeter(stim, value=rp.average_perimeter)

    elif prop == _VisProp.TOTAL_PERIMETER:
        return total_perimeter(stim, value=rp.total_perimeter)

    elif prop == _VisProp.AV_SURFACE_AREA:
        return average_surface_area(stim, value=rp.average_surface_area)

    elif prop == _VisProp.TOTAL_SURFACE_AREA:
        return total_surface_area(stim, value=rp.total_surface_area)

    elif prop == _VisProp.LOG_SIZE:
        return log_size(stim, value=rp.log_size)

    elif prop == _VisProp.LOG_SPACING:
        return log_spacing(stim, value=rp.log_spacing)

    elif prop == _VisProp.SPARSITY:
        return sparsity(stim, value=rp.sparsity)

    elif prop == _VisProp.FIELD_AREA:
        return field_area(stim, value=rp.field_area)

    elif prop == _VisProp.COVERAGE:
        return coverage(stim, value=rp.coverage)
    else:
        raise NotImplementedError(
            f"Not implemented for {str(prop)}"
        )


def property_difference(stim_pair: _NSNStimulusPair,
                        prop: _VisProp,
                        delta: _tp.Union[float, _np.float64],
                        adapt_stim: str = "both"):
    """Adapt visual property difference of NSNStimulusPair. Changes the property
    `prop` of the stimuli so that difference is equal to `delta`.

    There are three different adapt methods:
        `"a"`: change the properties of stimulus A only
        `"b"`: change the properties of stimulus B only
        `"both"`: change the stimulus A and B, each by 50% for the required change
    """
    p_a = stim_pair.stim_a.properties.get(prop)
    p_b = stim_pair.stim_b.properties.get(prop)
    rc = _np.float64(delta) - (p_a - p_b)  # required change

    if adapt_stim == "both":
        rc = rc / 2
        property_adapt(stim_pair.stim_a, prop, p_a + rc)
        property_adapt(stim_pair.stim_b, prop, p_b - rc)
    elif adapt_stim == "a":
        property_adapt(stim_pair.stim_a, prop, p_a + rc)
    elif adapt_stim == "b":
        property_adapt(stim_pair.stim_b, prop, p_b - rc)
    else:
        raise ValueError(f"Unknown adapt method {adapt_stim}. " +
                         "Must be either 'both', 'a' or 'b'")


def property_ratio(stim_pair: _NSNStimulusPair,
                   prop: _VisProp,
                   ratio: _tp.Union[float, _np.float64],
                   adapt_stim: str = "both"):
    """Adapt visual property ratio of NSNStimulusPair. Changes the property
    `prop` of the two stimuli so that ratio is equal to `delta`.

    There are three different adapt methods:
        `"a"`: change the properties of stimulus A only
        `"b"`: change the properties of stimulus B only
        `"both"`: change the stimulus A and B, each by 50% for the required change
    """

    pa = stim_pair.stim_a.properties.get(prop)
    pb = stim_pair.stim_b.properties.get(prop)

    rc = _np.float64(ratio) / (pa/pb)  # required change

    if adapt_stim == "both":
        rc = _np.sqrt(rc)
        property_adapt(stim_pair.stim_a, prop, pa * rc)
        property_adapt(stim_pair.stim_b, prop, pb / rc)
    elif adapt_stim == "a":
        property_adapt(stim_pair.stim_a, prop, pa * rc)
    elif adapt_stim == "b":
        property_adapt(stim_pair.stim_b, prop, pb / rc)
    else:
        raise ValueError(f"Unknown adapt method {property_adapt}. " +
                         "Must be either 'both', 'a' or 'b'")
