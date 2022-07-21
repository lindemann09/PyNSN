# calculates visual properties of a dot array/ dot cloud

from collections import OrderedDict
from enum import IntFlag, auto
from math import log2
import numpy as np

from .._lib.lib_typing import Any, Optional, NumPair

from .. import _arrays
from .convex_hull import ConvexHull, ConvexHullPositions


class VisualPropertyFlag(IntFlag):
    AV_DOT_DIAMETER = auto()
    AV_SURFACE_AREA = auto()
    AV_PERIMETER = auto()
    AV_RECT_SIZE = auto()

    TOTAL_SURFACE_AREA = auto()
    TOTAL_PERIMETER = auto()
    SPARSITY = auto()
    FIELD_AREA = auto()
    FIELD_AREA_POSITIONS = auto()
    COVERAGE = auto()

    LOG_SPACING = auto()
    LOG_SIZE = auto()

    NUMEROSITY = auto()

    def is_dependent_from(self, other_property: Any) -> bool:
        """returns true if both properties are not independent"""
        return (self.is_size_property() and other_property.is_size_property()) or \
               (self.is_space_property() and other_property.is_space_property())

    def is_size_property(self) -> bool :
        return self in (VisualPropertyFlag.LOG_SIZE,
                        VisualPropertyFlag.TOTAL_SURFACE_AREA,
                        VisualPropertyFlag.AV_DOT_DIAMETER,
                        VisualPropertyFlag.AV_SURFACE_AREA,
                        VisualPropertyFlag.AV_PERIMETER,
                        VisualPropertyFlag.TOTAL_PERIMETER)

    def is_space_property(self) -> bool:
        return self in (VisualPropertyFlag.LOG_SPACING,
                        VisualPropertyFlag.SPARSITY,
                        VisualPropertyFlag.FIELD_AREA)

    def label(self) -> str:
        labels = {
            VisualPropertyFlag.NUMEROSITY: "Numerosity",
            VisualPropertyFlag.LOG_SIZE: "Log Size",
            VisualPropertyFlag.TOTAL_SURFACE_AREA: "Total surface area",
            VisualPropertyFlag.AV_DOT_DIAMETER: "Average dot diameter",
            VisualPropertyFlag.AV_SURFACE_AREA: "Average surface area",
            VisualPropertyFlag.AV_PERIMETER: "Average perimeter",
            VisualPropertyFlag.TOTAL_PERIMETER: "Total perimeter",
            VisualPropertyFlag.AV_RECT_SIZE: "Average Rectangle Size",
            VisualPropertyFlag.LOG_SPACING: "Log Spacing",
            VisualPropertyFlag.SPARSITY: "Sparsity",
            VisualPropertyFlag.FIELD_AREA: "Field area",
            VisualPropertyFlag.FIELD_AREA_POSITIONS: "Field area positions",
            VisualPropertyFlag.COVERAGE: "Coverage"}
        return labels[self]


class ArrayProperties(object):

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
        if self._convex_hull is None:
            self._convex_hull = ConvexHull(self.oa)
        return self._convex_hull

    @property
    def convex_hull_positions(self) -> ConvexHullPositions:
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
    def average_rectangle_size(self) -> Optional[NumPair]:
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

    def get(self, property_flag: VisualPropertyFlag) -> Any:
        """returns a visual property"""

        assert isinstance(property_flag, VisualPropertyFlag)

        # Adapt
        if property_flag == VisualPropertyFlag.AV_DOT_DIAMETER:
            return self.average_dot_diameter

        elif property_flag == VisualPropertyFlag.AV_RECT_SIZE:
            return self.average_rectangle_size

        elif property_flag == VisualPropertyFlag.AV_PERIMETER:
            return self.average_perimeter

        elif property_flag == VisualPropertyFlag.TOTAL_PERIMETER:
            return self.total_perimeter

        elif property_flag == VisualPropertyFlag.AV_SURFACE_AREA:
            return self.average_surface_area

        elif property_flag == VisualPropertyFlag.TOTAL_SURFACE_AREA:
            return self.total_surface_area

        elif property_flag == VisualPropertyFlag.LOG_SIZE:
            return self.log_size

        elif property_flag == VisualPropertyFlag.LOG_SPACING:
            return self.log_spacing

        elif property_flag == VisualPropertyFlag.SPARSITY:
            return self.sparsity

        elif property_flag == VisualPropertyFlag.FIELD_AREA:
            return self.field_area

        elif property_flag == VisualPropertyFlag.FIELD_AREA_POSITIONS:
            return self.field_area_positions

        elif property_flag == VisualPropertyFlag.COVERAGE:
            return self.converage

        else:
            raise ValueError("{} is a unknown visual feature".format(property_flag))

    def as_dict(self) -> dict:
        """ordered dictionary with the most important visual properties"""
        rtn = [("Hash", self.oa.hash),
               ("Numerosity", self.numerosity),
               ("?", None),  # placeholder
               (VisualPropertyFlag.AV_PERIMETER.label(), self.average_perimeter),
               (VisualPropertyFlag.AV_SURFACE_AREA.label(), self.average_surface_area),
               (VisualPropertyFlag.TOTAL_PERIMETER.label(), self.total_perimeter),
               (VisualPropertyFlag.TOTAL_SURFACE_AREA.label(), self.total_surface_area),
               (VisualPropertyFlag.FIELD_AREA.label(), self.field_area),
               (VisualPropertyFlag.SPARSITY.label(), self.sparsity),
               (VisualPropertyFlag.COVERAGE.label(), self.converage),
               (VisualPropertyFlag.LOG_SIZE.label(), self.log_size),
               (VisualPropertyFlag.LOG_SPACING.label(), self.log_spacing)]

        if isinstance(self.oa, _arrays.DotArray):
            rtn[2] = (VisualPropertyFlag.AV_DOT_DIAMETER.label(), self.average_dot_diameter)
        elif isinstance(self.oa, _arrays.RectangleArray):
            rtn[2] = (VisualPropertyFlag.AV_RECT_SIZE.label(), self.average_rectangle_size)
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
