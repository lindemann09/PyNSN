# calculates visual properties of a nsn stimulus/ dot cloud

from collections import OrderedDict
from typing import Any, Optional, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ... import _stimulus, constants
from .._lib import geometry, rng
from .._lib.exceptions import NoSolutionError
from ._object_arrays.convex_hull import ConvexHull, ConvexHullPositions
from .. import VisualPropertyFlags


class ArrayProperties(object):
    """Visual properties fo an associated arrays

    Visual features of each nsn stimulus can be access and
    modified via an instance of this class
    """




    def fit_field_area(self, value: float, precision: Optional[float] = None) -> None:
        """changes the convex hull area to a desired size with certain precision

        uses scaling radial positions if field area has to be increased
        uses replacement of outer points (and later re-scaling)

        iterative method can takes some time.
        """

        if precision is None:
            precision = constants.DEFAULT_FIT_SPACING_PRECISION

        if self.field_area is None:
            return None  # not defined
        else:
            _match_field_area(self, value=value, precision=precision)

    def fit_coverage(
        self,
        value: float,
        precision: Optional[float] = None,
        fa2ta_ratio: Optional[float] = None,
    ) -> None:
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
            fa2ta_ratio = constants.DEFAULT_FIT_FA2TA_RATIO
        elif fa2ta_ratio < 0 or fa2ta_ratio > 1:
            fa2ta_ratio = 0.5
        if precision is None:
            precision = constants.DEFAULT_FIT_SPACING_PRECISION

        total_area_change100 = (value * self.field_area) - \
            self.total_surface_area
        d_change_total_area = total_area_change100 * (1 - fa2ta_ratio)
        if abs(d_change_total_area) > 0:
            self.fit_total_surface_area(
                self.total_surface_area + d_change_total_area)

        self.fit_field_area(self.total_surface_area /
                            value, precision=precision)





    def fit_log_spacing(self, value: float, precision: Optional[float] = None) -> None:
        """

        Parameters
        ----------
        value
        precision

        Returns
        -------

        """
        logfa = 0.5 * value + 0.5 * np.log2(self.numerosity)
        self.fit_field_area(value=2**logfa, precision=precision)

    def fit_log_size(self, value: float) -> None:
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        logtsa = 0.5 * value + 0.5 * np.log2(self.numerosity)
        self.fit_total_surface_area(2**logtsa)

    def fit_sparcity(self, value: float, precision=None) -> None:
        """

        Parameters
        ----------
        value
        precision

        Returns
        -------

        """
        return self.fit_field_area(value=value * self.numerosity, precision=precision)

    def fit(self, property_flag: VisualPropertyFlags, value: float) -> Any:
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
        if not isinstance(property_flag, VisualPropertyFlags):
            raise ValueError(
                "{} is not a visual feature.".format(property_flag))

        # Adapt
        if property_flag == VisualPropertyFlags.AV_DOT_DIAMETER:
            return self.fit_average_diameter(value=value)

        elif property_flag == VisualPropertyFlags.NUMEROSITY:
            return self.fit_numerosity(value=int(value))

        elif property_flag == VisualPropertyFlags.AV_PERIMETER:
            return self.fit_average_perimeter(value=value)

        elif property_flag == VisualPropertyFlags.TOTAL_PERIMETER:
            return self.fit_total_perimeter(value=value)

        elif property_flag == VisualPropertyFlags.AV_SURFACE_AREA:
            return self.fit_average_surface_area(value=value)

        elif property_flag == VisualPropertyFlags.TOTAL_SURFACE_AREA:
            return self.fit_total_surface_area(value=value)

        elif property_flag == VisualPropertyFlags.LOG_SIZE:
            return self.fit_log_size(value=value)

        elif property_flag == VisualPropertyFlags.LOG_SPACING:
            return self.fit_log_spacing(value=value)

        elif property_flag == VisualPropertyFlags.SPARSITY:
            return self.fit_sparcity(value=value)

        elif property_flag == VisualPropertyFlags.FIELD_AREA:
            return self.fit_field_area(value=value)

        elif property_flag == VisualPropertyFlags.COVERAGE:
            return self.fit_coverage(value=value)
        else:
            raise NotImplementedError(
                "Not implemented for {}".format(property_flag.label())
            )

    def scale_average_diameter(self, factor: float) -> None:
        if not isinstance(self._oa, DotArray):
            raise TypeError(
                "Scaling diameter is not possible for {}.".format(
                    type(self._oa).__name__
                )
            )
        if factor == 1:
            return
        return self.fit_average_diameter(self.average_dot_diameter * factor)

    def scale_average_rectangle_size(self, factor: float) -> None:
        if not isinstance(self._oa, RectangleArray):
            raise TypeError(
                "Scaling rectangle size is not possible for {}.".format(
                    type(self._oa).__name__
                )
            )
        if factor == 1:
            return
        return self.fit_average_rectangle_size(self.average_rectangle_size * factor)

    def scale_total_surface_area(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_total_surface_area(self.total_surface_area * factor)

    def scale_field_area(
        self, factor: float, precision: Optional[float] = None
    ) -> None:
        if factor == 1:
            return
        return self.fit_field_area(self.field_area * factor, precision=precision)

    def scale_coverage(
        self,
        factor: float,
        precision: Optional[float] = None,
        FA2TA_ratio: Optional[float] = None,
    ) -> None:
        if factor == 1:
            return
        return self.fit_coverage(
            self.coverage * factor, precision=precision, fa2ta_ratio=FA2TA_ratio
        )

    def scale_average_perimeter(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_average_perimeter(self.average_perimeter * factor)

    def scale_total_perimeter(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_total_perimeter(self.total_perimeter * factor)

    def scale_average_surface_area(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_average_surface_area(self.average_surface_area * factor)

    def scale_log_spacing(
        self, factor: float, precision: Optional[float] = None
    ) -> None:
        if factor == 1:
            return
        return self.fit_log_spacing(self.log_spacing * factor, precision=precision)

    def scale_log_size(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_log_size(self.log_size * factor)

    def scale_sparcity(self, factor: float, precision: Optional[float] = None) -> None:
        if factor == 1:
            return
        return self.fit_sparcity(self.sparsity * factor, precision=precision)

    def scale(self, feature: VisualPropertyFlags, factor: float) -> None:
        if factor == 1:
            return
        return self.fit(property_flag=feature, value=self.get(feature) * factor)


def _match_field_area(
    properties: ArrayProperties, value: float, precision: float
) -> None:
    """change the convex hull area to a desired size by scale the polar
    positions with certain precision

    iterative method can takes some time.

    Note: see doc string `field_area`
    """

    object_array = properties._oa
    current = properties.field_area

    if current is None:
        return  # not defined

    scale = 1  # find good scale
    step = 0.1
    if value < current:  # current too larger
        step *= -1

    # centered points
    old_center = properties.convex_hull.center
    object_array.xy = object_array.xy - old_center
    centered_polar = geometry.cartesian2polar(object_array.xy)

    # iteratively determine scale
    while abs(current - value) > precision:
        scale += step

        object_array.xy = geometry.polar2cartesian(centered_polar * [scale, 1])
        properties.reset()  # required at this point to recalc convex hull
        current = properties.field_area

        if (current < value and step < 0) or (current > value and step > 0):
            step *= -0.2  # change direction and finer grain

    object_array.xy = object_array.xy + old_center
    properties.reset()

    # TODO "visual test" (eye inspection) of fitting rect _arrays
