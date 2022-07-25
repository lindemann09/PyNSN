# calculates visual properties of a dot array/ dot cloud

from collections import OrderedDict
from math import log2
import numpy as np

from .. import _arrays
from .._lib.lib_typing import Any, Optional, NumPair, OptFloat, NDArray
from .._lib import rng
from .._lib import constants
from .._lib.constants import VisualPropertyFlags
from .._lib.exception import NoSolutionError
from .convex_hull import ConvexHull, ConvexHullPositions
from .tools import scale_field_area


class ArrayProperties(object):
    """Visual features of each object array can be access and
    modified via an instance of this class

    """


    def __init__(self, object_array: Any) -> None:
        # _lib or dot_cloud
        assert isinstance(object_array, (_arrays.DotArray, _arrays.RectangleArray,
                                         _arrays.PointArray))
        self.oa = object_array
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
            self._convex_hull = ConvexHull(self.oa)
        return self._convex_hull

    @property
    def convex_hull_positions(self) -> ConvexHullPositions:
        """TODO """
        if self._convex_hull_positions is None:
            self._convex_hull_positions = ConvexHullPositions(self.oa)
        return self._convex_hull_positions

    @property
    def average_dot_diameter(self) -> Optional[float]:
        if not isinstance(self.oa, _arrays.DotArray):
            return None
        elif self.numerosity == 0:
            return None
        else:
            return float(np.mean(self.oa.diameter))

    @property
    def average_rectangle_size(self) -> Optional[NDArray[float]]:
        if not isinstance(self.oa, _arrays.RectangleArray):
            return None
        elif self.numerosity == 0:
            return np.array([np.na, np.na])
        else:
            return np.mean(self.oa.sizes, axis=0)

    @property
    def total_surface_area(self) -> Optional[float]:
        return float(np.sum(self.oa.surface_areas))

    @property
    def average_surface_area(self) -> Optional[float]:
        if self.numerosity == 0:
            return None
        return float(np.mean(self.oa.surface_areas))

    @property
    def total_perimeter(self) -> Optional[float]:
        if self.numerosity == 0:
            return None
        return float(np.sum(self.oa.perimeter))

    @property
    def average_perimeter(self) -> Optional[float]:
        if self.numerosity == 0:
            return None
        return float(np.mean(self.oa.perimeter))

    @property
    def field_area_positions(self) -> Optional[float]:
        return self.convex_hull_positions.field_area

    @property
    def numerosity(self) -> int:
        return len(self.oa._xy)

    @property
    def converage(self) -> Optional[float]:
        """ percent coverage in the field area. It takes thus the object size
        into account. In contrast, the sparsity is only the ratio of field
        array and numerosity
        """
        try:
            return self.total_surface_area / self.field_area
        except ZeroDivisionError:
            return None

    @property
    def log_size(self) -> Optional[float]:
        try:
            return log2(self.total_surface_area) + \
                   log2(self.average_surface_area)
        except ValueError:
            return None

    @property
    def log_spacing(self) -> Optional[float]:
        try:
            return log2(self.field_area) + log2(self.sparsity)
        except ValueError:
            return None

    @property
    def sparsity(self) -> Optional[float]:
        try:
            return self.field_area / self.numerosity
        except ZeroDivisionError:
            return None

    @property
    def field_area(self) -> Optional[float]:
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
            return self.converage

        elif property_flag == VisualPropertyFlags.NUMEROSITY:
            return self.numerosity

        else:
            raise ValueError("{} is a unknown visual feature".format(property_flag))

    def as_dict(self) -> dict:
        """ordered dictionary with the most important visual properties"""
        rtn = [("Hash", self.oa.hash),
               ("Numerosity", self.numerosity),
               ("?", None),  # placeholder
               (VisualPropertyFlags.AV_PERIMETER.label(), self.average_perimeter),
               (VisualPropertyFlags.AV_SURFACE_AREA.label(), self.average_surface_area),
               (VisualPropertyFlags.TOTAL_PERIMETER.label(), self.total_perimeter),
               (VisualPropertyFlags.TOTAL_SURFACE_AREA.label(), self.total_surface_area),
               (VisualPropertyFlags.FIELD_AREA.label(), self.field_area),
               (VisualPropertyFlags.SPARSITY.label(), self.sparsity),
               (VisualPropertyFlags.COVERAGE.label(), self.converage),
               (VisualPropertyFlags.LOG_SIZE.label(), self.log_size),
               (VisualPropertyFlags.LOG_SPACING.label(), self.log_spacing)]

        if isinstance(self.oa, _arrays.DotArray):
            rtn[2] = (VisualPropertyFlags.AV_DOT_DIAMETER.label(), self.average_dot_diameter)
        elif isinstance(self.oa, _arrays.RectangleArray):
            rtn[2] = (VisualPropertyFlags.AV_RECT_SIZE.label(), self.average_rectangle_size)
        else:
            rtn.pop(2)
        return OrderedDict(rtn)

    def __str__(self) -> str:
        return self.as_text()

    def as_text(self, with_hash: bool = True,
                extended_format: bool = False,
                spacing_char: str = ".") -> str:
        if extended_format:
            rtn = None
            for k, v in self.as_dict().items():
                if rtn is None:
                    if with_hash:
                        rtn = "- {}: {}\n".format(k, v)
                    else:
                        rtn = ""
                else:
                    if rtn == "":
                        name = "- " + k
                    else:
                        name = "  " + k
                    try:
                        value = "{0:.2f}\n".format(v)  # try rounding
                    except:
                        value = "{}\n".format(v)

                    n_space = 14 - len(value)
                    if n_space < 2:
                        n_space = 2
                    rtn += name + (spacing_char * (24 - len(name))) + (" " * n_space) + value
        else:
            if with_hash:
                rtn = "ID: {} ".format(self.oa.hash)
            else:
                rtn = ""
            rtn += "N: {}, TSA: {}, ISA: {}, FA: {}, SPAR: {:.3f}, logSIZE: " \
                   "{:.2f}, logSPACE: {:.2f} COV: {:.2f}".format(
                self.numerosity,
                int(self.total_surface_area),
                int(self.average_surface_area),
                int(self.field_area),
                self.sparsity,
                self.log_size,
                self.log_spacing,
                self.converage)
        return rtn.rstrip()

    def fit_numerosity(self, value: int,
                       keep_convex_hull: bool = False) -> None:
        """

        """

        # make a copy for the deviant
        if value <= 0:
            self.oa.clear()
        else:
            # add or remove random dots
            change_numerosity = value - self.numerosity
            for _ in range(abs(change_numerosity)):
                if change_numerosity < 0:
                    # remove dots
                    if keep_convex_hull:
                        # find a random object that is not in convex hull
                        delete_id = None
                        ch = self.convex_hull.object_indices_unique
                        rnd_seq = list(range(0, self.numerosity))
                        rng.generator.shuffle(rnd_seq)
                        for x in rnd_seq:
                            if x not in ch:
                                delete_id = x
                                break
                        if delete_id is None:
                            raise NoSolutionError("Can't increase numeroisty, while keeping field area.")
                    else:
                        delete_id = rng.generator.integers(0, self.numerosity)

                    self.oa.delete(delete_id)

                else:
                    # add dot: copy a random dot
                    clone_id = rng.generator.integers(0, self.numerosity)
                    rnd_object = next(self.oa.iter_objects(clone_id))
                    try:
                        rnd_object = self.oa.get_free_position(
                            ref_object=rnd_object, allow_overlapping=False,
                            inside_convex_hull=keep_convex_hull)
                    except NoSolutionError:
                        # no free position
                        raise NoSolutionError("Can't increase numerosity. No free position found.")

                    self.oa.add([rnd_object])

    def fit_average_diameter(self, value: float) -> None:
        # changes diameter
        if not isinstance(self.oa, _arrays.DotArray):
            raise TypeError("Adapting diameter is not possible for {}.".format(
                type(self.oa).__name__))
        scale = value / self.average_dot_diameter
        self.oa._diameter = self.oa.diameter * scale
        self.reset()

    def fit_average_rectangle_size(self, value: NumPair) -> None:
        # changes diameter
        if not isinstance(self.oa, _arrays.RectangleArray):
            raise TypeError("Adapting rectangle size is not possible for {}.".format(
                type(self.oa).__name__))
        try:
            width, height = value
        except TypeError:
            raise TypeError("Value ({}) has to tuple of 2 numerical (width, height).".format(
                value))

        scale = np.array([width / self.average_rectangle_size[0],
                          height / self.average_rectangle_size[1]])
        self.oa._sizes = self.oa._sizes * scale
        self.reset()

    def fit_total_surface_area(self, value: float) -> None:
        """
        
        Parameters
        ----------
        value

        Returns
        -------

        """
        a_scale = value / self.total_surface_area
        if isinstance(self.oa, _arrays.DotArray):
            self.oa._diameter = np.sqrt(
                self.oa.surface_areas * a_scale) * 2 / np.sqrt(
                np.pi)  # d=sqrt(4a/pi) = sqrt(a)*2/sqrt(pi)
        else:  # rect
            self.oa._sizes = self.oa._sizes * np.sqrt(a_scale)

        self.reset()

    def fit_field_area(self, value: float,
                       precision: OptFloat = None) -> None:
        """
        
        Parameters
        ----------
        value
        precision

        Returns
        -------

        """
        
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
            scale_field_area(self.oa, value=value, precision=precision)

    def fit_coverage(self, value: float,
                 precision: OptFloat = None,
                 FA2TA_ratio: OptFloat = None) -> None:
        """
        
        Parameters
        ----------
        value
        precision
        FA2TA_ratio

        Returns
        -------

        """

        # FIXME check drifting outwards if extra space is small and adapt_FA2TA_ratio=1
        # when to realign, realignment changes field_area!
        """this function changes the area and remixes to get a desired density
        precision in percent between 1 < 0

        ratio_area_convex_hull_adaptation:
            ratio of adaptation via area or via convex_hull (between 0 and 1)

        """
        
        print("WARNING: _adapt_coverage is a experimental ")
        # dens = convex_hull_area / total_surface_area
        if FA2TA_ratio is None:
            FA2TA_ratio = constants.DEFAULT_FIT_FA2TA_RATIO
        elif FA2TA_ratio < 0 or FA2TA_ratio > 1:
            FA2TA_ratio = 0.5
        if precision is None:
            precision = constants.DEFAULT_FIT_SPACING_PRECISION

        total_area_change100 = (value * self.field_area) -  self.total_surface_area
        d_change_total_area = total_area_change100 * (1 - FA2TA_ratio)
        if abs(d_change_total_area) > 0:
            self.fit_total_surface_area(self.total_surface_area + d_change_total_area)

        self.fit_field_area(self.total_surface_area / value, precision=precision)

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
        if isinstance(self.oa, _arrays.DotArray):
            self.fit_average_diameter(value / (self.numerosity * np.pi))

        elif isinstance(self.oa, _arrays.RectangleArray):
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
                        precision: OptFloat = None) -> None:
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
        self.fit_total_surface_area( 2 ** logtsa)

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
            raise ValueError("{} is not a visual feature.".format(property_flag))

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
        if not isinstance(self.oa, _arrays.DotArray):
            raise TypeError("Scaling diameter is not possible for {}.".format(
                type(self.oa).__name__))
        if factor == 1:
            return
        return self.fit_average_diameter(self.average_dot_diameter * factor)

    def scale_average_rectangle_size(self, factor: float) -> None:
        if not isinstance(self.oa, _arrays.RectangleArray):
            raise TypeError("Scaling rectangle size is not possible for {}.".format(
                type(self.oa).__name__))
        if factor == 1:
            return
        return self.fit_average_rectangle_size(self.average_rectangle_size * factor)

    def scale_total_surface_area(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_total_surface_area(self.total_surface_area * factor)

    def scale_field_area(self, factor: float,
                         precision: OptFloat = None) -> None:
        if factor == 1:
            return
        return self.fit_field_area(self.field_area * factor, precision=precision)

    def scale_coverage(self, factor: float,
                       precision: OptFloat = None,
                       FA2TA_ratio: OptFloat = None) -> None:
        if factor == 1:
            return
        return self.fit_coverage(self.converage * factor,
                                 precision=precision,
                                 FA2TA_ratio=FA2TA_ratio)

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
                          precision: OptFloat = None) -> None:
        if factor == 1:
            return
        return self.fit_log_spacing(self.log_spacing * factor, precision=precision)

    def scale_log_size(self, factor: float) -> None:
        if factor == 1:
            return
        return self.fit_log_size(self.log_size * factor)

    def scale_sparcity(self, factor: float,
                       precision: OptFloat = None) -> None:
        if factor == 1:
            return
        return self.fit_sparcity(self.sparsity * factor, precision=precision)

    def scale(self, feature: VisualPropertyFlags, factor: float) -> None:
        if factor == 1:
            return
        return self.fit(property_flag=feature,
                                        value=self.get(feature) * factor)

    # TODO "visual test" (eye inspection) of fitting rect _arrays
