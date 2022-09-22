# calculates visual properties of a dot array/ dot cloud

from collections import OrderedDict
from typing import Any, Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .. import _arrays, constants
from .._lib import geometry, rng
from ..constants import VisualPropertyFlags
from .._lib.exception import NoSolutionError
from .convex_hull import ConvexHull, ConvexHullPositions


class ArrayProperties(object):
    """Visual properties fo an associated arrays

    Visual features of each object array can be access and
    modified via an instance of this class
    """

    def __init__(self, object_array: Any) -> None:
        # _lib or dot_cloud
        assert isinstance(object_array, (_arrays.DotArray,
                                         _arrays.RectangleArray))
        self._oa = object_array
        self._convex_hull = None
        self._convex_hull_positions = None

    def reset(self) -> None:
        """reset to enforce recalculation of certain parameter """
        self._convex_hull = None
        self._convex_hull_positions = None

    @property
    def convex_hull(self) -> ConvexHull:
        """TODO """
        if self._convex_hull is None:
            self._convex_hull = ConvexHull(self._oa)
        return self._convex_hull

    @property
    def convex_hull_positions(self) -> ConvexHullPositions:
        """TODO """
        if self._convex_hull_positions is None:
            self._convex_hull_positions = ConvexHullPositions(self._oa)
        return self._convex_hull_positions

    @property
    def average_dot_diameter(self) -> float:
        if not isinstance(self._oa, _arrays.DotArray) or self.numerosity == 0:
            return np.nan
        else:
            return float(np.mean(self._oa.diameter))

    @property
    def average_rectangle_size(self) -> NDArray:
        if not isinstance(self._oa, _arrays.RectangleArray) \
                or self.numerosity == 0:
            return np.array([np.nan, np.nan])
        else:
            return np.mean(self._oa.sizes, axis=0)

    @property
    def total_surface_area(self) -> float:
        return float(np.sum(self._oa.surface_areas))

    @property
    def average_surface_area(self) -> float:
        if self.numerosity == 0:
            return 0
        return float(np.mean(self._oa.surface_areas))

    @property
    def total_perimeter(self) -> float:
        if self.numerosity == 0:
            return 0
        return float(np.sum(self._oa.perimeter))

    @property
    def average_perimeter(self) -> float:
        if self.numerosity == 0:
            return np.nan
        return float(np.mean(self._oa.perimeter))

    @property
    def field_area_positions(self) -> float:
        return self.convex_hull_positions.field_area

    @property
    def numerosity(self) -> int:
        return len(self._oa.xy)

    @property
    def coverage(self) -> float:
        """ percent coverage in the field area. It takes thus the object size
        into account. In contrast, the sparsity is only the ratio of field
        array and numerosity
        """
        try:
            return self.total_surface_area / self.field_area
        except ZeroDivisionError:
            return np.nan

    @property
    def log_size(self) -> float:
        try:
            return np.log2(self.total_surface_area) + \
                np.log2(self.average_surface_area)
        except ValueError:
            return np.nan

    @property
    def log_spacing(self) -> float:
        try:
            return np.log2(self.field_area) + np.log2(self.sparsity)
        except ValueError:
            return np.nan

    @property
    def sparsity(self) -> float:
        try:
            return self.field_area / self.numerosity
        except ZeroDivisionError:
            return np.nan

    @property
    def field_area(self) -> float:
        return self.convex_hull.field_area

    def get(self, property_flag: VisualPropertyFlags) -> Any:
        """returns a visual property"""

        assert isinstance(property_flag, VisualPropertyFlags)

        # Adapt
        if property_flag == VisualPropertyFlags.AV_DOT_DIAMETER:
            return self.average_dot_diameter

        elif property_flag == VisualPropertyFlags.AV_RECT_SIZE:
            return self.average_rectangle_size

        elif property_flag == VisualPropertyFlags.AV_PERIMETER:
            return self.average_perimeter

        elif property_flag == VisualPropertyFlags.TOTAL_PERIMETER:
            return self.total_perimeter

        elif property_flag == VisualPropertyFlags.AV_SURFACE_AREA:
            return self.average_surface_area

        elif property_flag == VisualPropertyFlags.TOTAL_SURFACE_AREA:
            return self.total_surface_area

        elif property_flag == VisualPropertyFlags.LOG_SIZE:
            return self.log_size

        elif property_flag == VisualPropertyFlags.LOG_SPACING:
            return self.log_spacing

        elif property_flag == VisualPropertyFlags.SPARSITY:
            return self.sparsity

        elif property_flag == VisualPropertyFlags.FIELD_AREA:
            return self.field_area

        elif property_flag == VisualPropertyFlags.FIELD_AREA_POSITIONS:
            return self.field_area_positions

        elif property_flag == VisualPropertyFlags.COVERAGE:
            return self.coverage

        elif property_flag == VisualPropertyFlags.NUMEROSITY:
            return self.numerosity

        else:
            raise ValueError("f{property_flag} is a unknown visual feature")

    def to_dict(self) -> dict:
        """Dictionary with the visual properties"""
        rtn = [("Hash", self._oa.hash),
               ("Numerosity", self.numerosity),
               ("?", None),  # placeholder
               (VisualPropertyFlags.AV_PERIMETER.label(), self.average_perimeter),
               (VisualPropertyFlags.AV_SURFACE_AREA.label(),
                self.average_surface_area),
               (VisualPropertyFlags.TOTAL_PERIMETER.label(), self.total_perimeter),
               (VisualPropertyFlags.TOTAL_SURFACE_AREA.label(),
                self.total_surface_area),
               (VisualPropertyFlags.FIELD_AREA.label(), self.field_area),
               (VisualPropertyFlags.SPARSITY.label(), self.sparsity),
               (VisualPropertyFlags.COVERAGE.label(), self.coverage),
               (VisualPropertyFlags.LOG_SIZE.label(), self.log_size),
               (VisualPropertyFlags.LOG_SPACING.label(), self.log_spacing)]

        if isinstance(self._oa, _arrays.DotArray):
            rtn[2] = (VisualPropertyFlags.AV_DOT_DIAMETER.label(),
                      self.average_dot_diameter)
        elif isinstance(self._oa, _arrays.RectangleArray):
            rtn[2] = (VisualPropertyFlags.AV_RECT_SIZE.label(),
                      self.average_rectangle_size.tolist())
        else:
            rtn.pop(2)
        return OrderedDict(rtn)

    def __str__(self) -> str:
        return self.as_text(extended_format=True)

    def as_text(self, with_hash: bool = True,
                extended_format: bool = False,
                spacing_char: str = " ") -> str:
        if extended_format:
            rtn = None
            for k, v in self.to_dict().items():
                if rtn is None:
                    if with_hash:
                        rtn = f"- {k}: {v}\n"
                    else:
                        rtn = ""
                else:
                    if rtn == "":
                        name = "- " + k
                    else:
                        name = "  " + k
                    try:
                        value = f"{v:.2f}\n"  # try rounding
                    except (ValueError, TypeError):
                        value = f"{v}\n"

                    n_space = 14 - len(value)
                    if n_space < 2:
                        n_space = 2
                    rtn += name + (spacing_char * (24 - len(name))
                                   ) + (" " * n_space) + value
            if rtn is None:
                rtn = ""
        else:
            if with_hash:
                rtn = "ID: {} ".format(self._oa.hash)
            else:
                rtn = ""
            rtn += f"N: {self.numerosity}, " \
                + f"TSA: {int(self.total_surface_area)}, " \
                + f"ISA: {int(self.average_surface_area)}, "\
                + f"FA: {int(self.field_area)}, "\
                + f"SPAR: {self.sparsity:.3f}, "\
                + f"logSIZE: {self.log_size:.2f}, "\
                + f"logSPACE: {self.log_spacing:.2f}, "\
                + f"COV: {self.coverage:.2f}"

        return rtn.rstrip()

    def fit_numerosity(self, value: int,
                       keep_convex_hull: bool = False) -> None:
        """

        """

        # make a copy for the deviant
        if value <= 0:
            self._oa.clear()
        else:
            # add or remove random dots
            change_numerosity = value - self.numerosity
            for _ in range(abs(change_numerosity)):
                if change_numerosity < 0:
                    # remove dots
                    if keep_convex_hull:
                        # find a random object that is not in convex hull
                        delete_id = None
                        ch = self.convex_hull.object_indices
                        rnd_seq = list(range(0, self.numerosity))
                        rng.generator.shuffle(rnd_seq)
                        for x in rnd_seq:
                            if x not in ch:
                                delete_id = x
                                break
                        if delete_id is None:
                            raise NoSolutionError(
                                "Can't increase numeroisty, while keeping field area.")
                    else:
                        delete_id = rng.generator.integers(0, self.numerosity)

                    self._oa.delete(delete_id)

                else:
                    # add object: copy a random dot
                    clone_id = rng.generator.integers(0, self.numerosity)
                    rnd_object = next(self._oa.iter_objects(clone_id))
                    try:
                        rnd_object = self._oa.get_free_position(
                            ref_object=rnd_object, allow_overlapping=False,
                            inside_convex_hull=keep_convex_hull)
                    except NoSolutionError as err:
                        # no free position
                        raise NoSolutionError(
                            "Can't increase numerosity. No free position found.") from err

                    self._oa.add([rnd_object])

    def fit_average_diameter(self, value: float) -> None:
        """Set average diameter.

        Args:
            value: diameter

        Raises:
            TypeError: if associated array is not a DotArray
        """
        # changes diameter
        if not isinstance(self._oa, _arrays.DotArray):
            raise TypeError("Adapting diameter is not possible "
                            + f"for {type(self._oa).__name__}.")
        scale = value / self.average_dot_diameter
        self._oa._np_shapes.diameter = self._oa.diameter * scale  # pylint: disable=W0212
        self.reset()

    def fit_average_rectangle_size(self, value: ArrayLike) -> None:
        """Set average rectangle size.

        Args:
            value:  (width, height)

        Raises:
            TypeError: if associated array is not a RectangleArray or values is not a tuple of two numerical values
        """
        # changes diameter
        if not isinstance(self._oa, _arrays.RectangleArray):
            raise RuntimeError("Adapting rectangle size is not possible for "
                               + f"{type(self._oa).__name__}.")
        new_size = np.asarray(value)
        if new_size.shape != (2,):
            raise ValueError(f"Value ({value}) has to tuple of 2 numerical "
                             + "(width, height).")

        av_size = self.average_rectangle_size
        if not np.all(av_size > 0):
            raise RuntimeError(
                "Numerosity, width or hight is zero or not defined.")
        scale = np.divide(new_size, av_size)
        self._oa._np_shapes.sizes = self._oa._np_shapes.sizes * \
            scale  # pylint: disable=W0212
        self.reset()

    def fit_total_surface_area(self, value: float) -> None:
        """Set surface area.

        Resize all object to fit a specific surface area

        Args:
            value: surface area
        """
        a_scale = value / self.total_surface_area
        # pylint: disable=W0212
        if isinstance(self._oa, _arrays.DotArray):
            self._oa._np_shapes.diameter = np.sqrt(self._oa.surface_areas * a_scale) \
                * 2 / np.sqrt(np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
        elif isinstance(self._oa, _arrays.RectangleArray):
            # rect
            self._oa._np_shapes.sizes = self._oa._np_shapes.sizes * \
                np.sqrt(a_scale)

        self.reset()

    def fit_field_area(self, value: float,
                       precision: Optional[float] = None) -> None:
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
            _match_field_area(self._oa, value=value, precision=precision)

    def fit_coverage(self, value: float,
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

    def fit_average_perimeter(self, value: float) -> None:
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        self.fit_total_perimeter(value * self.numerosity)

    def fit_total_perimeter(self, value: float) -> None:
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        if isinstance(self._oa, _arrays.DotArray):
            self.fit_average_diameter(value / (self.numerosity * np.pi))

        elif isinstance(self._oa, _arrays.RectangleArray):
            new_size = self.average_rectangle_size * value / self.total_perimeter

            self.fit_average_rectangle_size(new_size)

    def fit_average_surface_area(self, value: float) -> None:
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        self.fit_total_surface_area(self.numerosity * value)

    def fit_log_spacing(self, value: float,
                        precision: Optional[float] = None) -> None:
        """

        Parameters
        ----------
        value
        precision

        Returns
        -------

        """
        logfa = 0.5 * value + 0.5 * np.log2(self.numerosity)
        self.fit_field_area(value=2 ** logfa, precision=precision)

    def fit_log_size(self, value: float) -> None:
        """

        Parameters
        ----------
        value

        Returns
        -------

        """
        logtsa = 0.5 * value + 0.5 * np.log2(self.numerosity)
        self.fit_total_surface_area(2 ** logtsa)

    def fit_sparcity(self, value: float,
                     precision=None) -> None:
        """

        Parameters
        ----------
        value
        precision

        Returns
        -------

        """
        return self.fit_field_area(value=value * self.numerosity,
                                   precision=precision)

    def fit(self, property_flag: VisualPropertyFlags,
            value: float) -> Any:
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
            raise NotImplementedError("Not implemented for {}".format(
                property_flag.label()))

    def scale_average_diameter(self, factor: float) -> None:
        if not isinstance(self._oa, _arrays.DotArray):
            raise TypeError("Scaling diameter is not possible for {}.".format(
                type(self._oa).__name__))
        if factor == 1:
            return
        return self.fit_average_diameter(self.average_dot_diameter * factor)

    def scale_average_rectangle_size(self, factor: float) -> None:
        if not isinstance(self._oa, _arrays.RectangleArray):
            raise TypeError("Scaling rectangle size is not possible for {}.".format(
                type(self._oa).__name__))
        if factor == 1:
            return
        return self.fit_average_rectangle_size(self.average_rectangle_size * factor)

    def scale_total_surface_area(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_total_surface_area(self.total_surface_area * factor)

    def scale_field_area(self, factor: float,
                         precision: Optional[float] = None) -> None:
        if factor == 1:
            return
        return self.fit_field_area(self.field_area * factor, precision=precision)

    def scale_coverage(self, factor: float,
                       precision: Optional[float] = None,
                       FA2TA_ratio: Optional[float] = None) -> None:
        if factor == 1:
            return
        return self.fit_coverage(self.coverage * factor,
                                 precision=precision,
                                 fa2ta_ratio=FA2TA_ratio)

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

    def scale_log_spacing(self, factor: float,
                          precision: Optional[float] = None) -> None:
        if factor == 1:
            return
        return self.fit_log_spacing(self.log_spacing * factor, precision=precision)

    def scale_log_size(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_log_size(self.log_size * factor)

    def scale_sparcity(self, factor: float,
                       precision: Optional[float] = None) -> None:
        if factor == 1:
            return
        return self.fit_sparcity(self.sparsity * factor, precision=precision)

    def scale(self, feature: VisualPropertyFlags, factor: float) -> None:
        if factor == 1:
            return
        return self.fit(property_flag=feature,
                        value=self.get(feature) * factor)


def _match_field_area(object_array,
                      value: float,
                      precision: float) -> None:
    """change the convex hull area to a desired size by scale the polar
    positions with certain precision

    iterative method can takes some time.

    Note: see doc string `field_area`
    """
    current = object_array.properties.field_area

    if current is None:
        return  # not defined

    scale = 1  # find good scale
    step = 0.1
    if value < current:  # current too larger
        step *= -1

    # centered points
    old_center = object_array.get_center_of_field_area()
    object_array._np_shapes.xy = object_array.xy - old_center
    centered_polar = geometry.cartesian2polar(object_array.xy)

    # iteratively determine scale
    while abs(current - value) > precision:
        scale += step

        object_array._np_shapes.xy = geometry.polar2cartesian(
            centered_polar * [scale, 1])
        object_array.properties.reset()  # required at this point to recalc convex hull
        current = object_array.properties.field_area

        if (current < value and step < 0) or \
                (current > value and step > 0):
            step *= -0.2  # change direction and finer grain

    object_array._np_shapes.xy = object_array.xy + old_center
    object_array.properties.reset()

    # TODO "visual test" (eye inspection) of fitting rect _arrays
