"""Fitting module"""
# pylint: disable=W0212
from typing import Any, Optional, Union

import numpy as np
import shapely

from pynsn.errors import NoSolutionError

from . import defaults
from ._misc import cartesian2polar, polar2cartesian
from ._stimulus import NSNStimulus, VisProp
from .random import generator

# FIXME module not yer tested


def total_surface_area(stim: NSNStimulus, value: Union[float, np.float_]) -> None:
    """Set surface area.

    Resize all object to fit a specific surface area

    Args:
        value: surface area
    """
    stim.scale(value / stim.properties.total_surface_area)


def total_perimeter(stim: NSNStimulus, value: Union[float, np.float_]) -> None:
    """fit the total parameter of the stimulus"""
    stim.scale(value / stim.properties.total_perimeter)


def average_perimeter(stim: NSNStimulus, value: Union[float, np.float_]) -> None:
    """fit the average parameter of the stimulus"""
    total_perimeter(stim, value * stim.n_objects)


def average_surface_area(stim: NSNStimulus, value: Union[float, np.float_]) -> None:
    """fits the average surface area of the stimulus"""
    total_surface_area(stim, stim.n_objects * value)


def numerosity(stim: NSNStimulus, value: int,
               keep_convex_hull: bool = False,
               max_iterations: Optional[int] = None) -> None:
    """
    fitting the numerosity
    """

    # make a copy for the deviant
    if value <= 0:
        stim.clear()
    else:
        # add or remove random dots
        change_numerosity = value - stim.n_objects
        if keep_convex_hull and change_numerosity < 0:
            # find objects touching the convex hull (ch_shapes)
            ring = shapely.get_exterior_ring(stim.convex_hull.polygon)
            ch_shapes = np.flatnonzero(shapely.intersects(
                stim.polygons, ring))
        else:
            ch_shapes = None

        for _ in range(abs(change_numerosity)):
            if change_numerosity < 0:
                # remove dots
                if ch_shapes is not None:
                    # find a random object that is not in convex hull
                    rnd_seq = np.arange(stim.n_objects)
                    generator.shuffle(rnd_seq)
                    delete_id = None
                    for x in rnd_seq:
                        if x not in ch_shapes:
                            delete_id = x
                            break
                    if delete_id is None:
                        raise NoSolutionError(
                            "Can't increase numerosity, while keeping field area.")
                else:
                    delete_id = generator.integers(0, stim.n_objects)

                stim.delete(delete_id)

            else:
                # add dot: copy a random dot
                clone_id = generator.integers(0, stim.n_objects)
                rnd_object = stim.shapes[clone_id]
                try:
                    rnd_object = stim.random_free_position(
                        shape=rnd_object, ignore_overlaps=False,
                        inside_convex_hull=keep_convex_hull,
                        max_iterations=max_iterations)
                except NoSolutionError as err:
                    # no free position
                    raise NoSolutionError(
                        "Can't increase numerosity. No free position found.") from err

                stim.add(rnd_object)


def field_area(stim: NSNStimulus, value: Union[float, np.float_],
               precision: Optional[Union[float, np.float_]] = None) -> None:
    """changes the convex hull area to a desired size with certain precision

    uses scaling radial positions if field area has to be increased
    uses replacement of outer points (and later re-scaling)

    iterative method can takes some time.
    """

    if precision is None or np.isnan(precision):
        precision = defaults.FIT_SPACING_PRECISION

    current = stim.convex_hull.area
    if stim.n_objects < 3 or current == 0:
        return None  # not defined

    scale = 1  # find good scale
    step = 0.1
    if value < current:  # current too larger
        step *= -1

    # centred points
    old_center = stim.convex_hull.centroid
    stim.set_xy(stim._xy - old_center)
    centered_polar = cartesian2polar(stim._xy)

    # iteratively determine scale
    while abs(current - value) > precision:
        scale += step

        stim.set_xy(polar2cartesian(centered_polar * [scale, 1]))
        current = stim.convex_hull.area

        if (current < value and step < 0) or (current > value and step > 0):
            step *= -0.2  # change direction and finer grain

    stim.set_sizes(stim._xy + old_center)  # move back


def coverage(stim: NSNStimulus, value: Union[float, np.float_],
             precision: Optional[float] = None,
             fa2ta_ratio: Optional[float] = None) -> None:
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
        fa2ta_ratio = defaults.FIT_FA2TA_RATIO
    elif fa2ta_ratio < 0 or fa2ta_ratio > 1:
        fa2ta_ratio = 0.5
    if precision is None:
        precision = defaults.FIT_SPACING_PRECISION

    ta = stim.properties.total_surface_area  # total area
    ta_change100 = (value * stim.properties.field_area) - ta
    d_ta_change = ta_change100 * (1 - fa2ta_ratio)
    if abs(d_ta_change) > 0:
        total_surface_area(stim, ta + d_ta_change)

    field_area(stim,
               value=stim.properties.total_surface_area / value,
               precision=precision)


def log_spacing(stim: NSNStimulus,
                value: Union[float, np.float_],
                precision: Optional[float] = None) -> None:
    """

    Parameters
    ----------
    value
    precision

    Returns
    -------

    """
    log_fa = 0.5 * value + 0.5 * np.log2(stim.n_objects)
    field_area(stim, value=2**log_fa, precision=precision)


def log_size(stim: NSNStimulus, value: Union[float, np.float_]) -> None:
    """

    Parameters
    ----------
    value

    Returns
    -------

    """
    log_tsa = 0.5 * value + 0.5 * np.log2(stim.n_objects)
    total_surface_area(stim, value=2**log_tsa)


def sparsity(stim: NSNStimulus,
             value: Union[float, np.float_], precision=None) -> None:
    """

    Parameters
    ----------
    value
    precision

    Returns
    -------

    """
    field_area(stim, value=value * stim.n_objects, precision=precision)


def fit(stim: NSNStimulus, property_flag: VisProp, value: float) -> Any:
    """
    adapt_properties: continuous property or list of continuous properties
    several properties to be adapted
    if adapt nsn stimulus is specified, array will be adapt to adapt_dot_array, otherwise
    the values defined in adapt_properties is used.
    some adapting requires realignement to avoid overlaps. However,
    realigment might result in a different field area. Thus, realign after
    adapting for  Size parameter and realign before adapting space
    parameter.

    """

    # type check
    if not isinstance(property_flag, VisProp):
        raise ValueError(
            f"{property_flag} is not a visual feature.")

    # Adapt
    if property_flag == VisProp.NUMEROSITY:
        return numerosity(stim, value=int(value))

    elif property_flag == VisProp.AV_PERIMETER:
        return average_perimeter(stim, value=value)

    elif property_flag == VisProp.TOTAL_PERIMETER:
        return total_perimeter(stim, value=value)

    elif property_flag == VisProp.AV_SURFACE_AREA:
        return average_surface_area(stim, value=value)

    elif property_flag == VisProp.TOTAL_SURFACE_AREA:
        return total_surface_area(stim, value=value)

    elif property_flag == VisProp.LOG_SIZE:
        return log_size(stim, value=value)

    elif property_flag == VisProp.LOG_SPACING:
        return log_spacing(stim, value=value)

    elif property_flag == VisProp.SPARSITY:
        return sparsity(stim, value=value)

    elif property_flag == VisProp.FIELD_AREA:
        return field_area(stim, value=value)

    elif property_flag == VisProp.COVERAGE:
        return coverage(stim, value=value)
    else:
        raise NotImplementedError(
            f"Not implemented for {property_flag.label()}"
        )
