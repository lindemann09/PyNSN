# calculates visual features of a dot array/ dot cloud

from collections import OrderedDict
from enum import IntFlag, auto

import numpy as np
from . import misc
from . import arrays
from .convex_hull import ConvexHull, ConvexHullPositions

class VisualFeature(IntFlag):

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

    def is_dependent_from(self, featureB):
        """returns true if both features are not independent"""
        return (self.is_size_feature() and featureB.is_size_feature()) or \
               (self.is_space_feature() and featureB.is_space_feature())

    def is_size_feature(self):
        return self in (VisualFeature.LOG_SIZE,
                        VisualFeature.TOTAL_SURFACE_AREA,
                        VisualFeature.AV_DOT_DIAMETER,
                        VisualFeature.AV_SURFACE_AREA,
                        VisualFeature.AV_PERIMETER,
                        VisualFeature.TOTAL_PERIMETER)

    def is_space_feature(self):
        return self in (VisualFeature.LOG_SPACING,
                        VisualFeature.SPARSITY,
                        VisualFeature.FIELD_AREA)

    def label(self):
        labels = {
            VisualFeature.LOG_SIZE: "Log Size",
            VisualFeature.TOTAL_SURFACE_AREA: "Total surface area",
            VisualFeature.AV_DOT_DIAMETER: "Average dot diameter",
            VisualFeature.AV_SURFACE_AREA: "Average surface area",
            VisualFeature.AV_PERIMETER: "Average perimeter",
            VisualFeature.TOTAL_PERIMETER: "Total perimeter",
            VisualFeature.AV_RECT_SIZE: "Average Rectangle Size",
            VisualFeature.LOG_SPACING: "Log Spacing",
            VisualFeature.SPARSITY: "Sparsity",
            VisualFeature.FIELD_AREA: "Field area",
            VisualFeature.COVERAGE: "Coverage"}
        return labels[self]


class ArrayFeatures(object):

    def __init__(self, object_array):
        # _lib or dot_cloud
        assert isinstance(object_array, arrays.GenericObjectArray)
        self.oa = object_array
        self._convex_hull = None
        self._convex_hull_positions = None

    def reset(self):
        """reset to enforce recalculation of certain parameter """
        self._convex_hull = None
        self._convex_hull_positions = None

    @property
    def convex_hull(self):
        if self._convex_hull is None:
            self._convex_hull = ConvexHull(self.oa)
        return self._convex_hull

    @property
    def convex_hull_positions(self):
        if self._convex_hull_positions is None:
            self._convex_hull_positions = ConvexHullPositions(self.oa)
        return self._convex_hull_positions

    @property
    def average_dot_diameter(self):
        if not isinstance(self.oa, arrays.DotArray):
            return None
        elif self.numerosity == 0:
            return np.nan
        else:
            return np.mean(self.oa.diameters)

    @property
    def average_rectangle_size(self):
        if not isinstance(self.oa, arrays.RectangleArray):
            return None
        elif self.numerosity == 0:
            return np.array([np.nan, np.nan])
        else:
            return np.mean(self.oa.sizes, axis=0)

    @property
    def total_surface_area(self):
        return np.sum(self.oa.surface_areas)

    @property
    def average_surface_area(self):
        if self.numerosity == 0:
            return np.nan
        return np.mean(self.oa.surface_areas)

    @property
    def total_perimeter(self):
        return np.sum(self.oa.perimeter)

    @property
    def average_perimeter(self):
        if self.numerosity == 0:
            return np.nan
        return np.mean(self.oa.perimeter)

    @property
    def field_area_positions(self):
        return self.convex_hull_positions.field_area

    @property
    def numerosity(self):
        return len(self.oa._xy)

    @property
    def converage(self):
        """ percent coverage in the field area. It takes thus the object size
        into account. In contrast, the sparsity is only the ratio of field
        array and numerosity
        """
        try:
            return self.total_surface_area / self.field_area
        except ZeroDivisionError:
            return np.nan

    @property
    def log_size(self):
        try:
            return misc.log2(self.total_surface_area) + misc.log2(
                    self.average_surface_area)
        except ValueError:
            return np.nan


    @property
    def log_spacing(self):
        try:
            return misc.log2(self.field_area) + misc.log2(self.sparsity)
        except ValueError:
            return np.nan

    @property
    def sparsity(self):
        try:
            return self.field_area / self.numerosity
        except ZeroDivisionError:
            return np.nan

    @property
    def field_area(self):
        return self.convex_hull.field_area

    def get(self, feature):
        """returns a feature"""

        assert isinstance(feature, VisualFeature)

       # Adapt
        if feature == VisualFeature.AV_DOT_DIAMETER:
            return self.average_dot_diameter

        elif feature == VisualFeature.AV_PERIMETER:
            return self.average_perimeter

        elif feature == VisualFeature.TOTAL_PERIMETER:
            return self.total_perimeter

        elif feature == VisualFeature.AV_SURFACE_AREA:
            return self.average_surface_area

        elif feature == VisualFeature.TOTAL_SURFACE_AREA:
            return self.total_surface_area

        elif feature == VisualFeature.LOG_SIZE:
            return self.log_size

        elif feature == VisualFeature.LOG_SPACING:
            return self.log_spacing

        elif feature == VisualFeature.SPARSITY:
            return self.sparsity

        elif feature == VisualFeature.FIELD_AREA:
            return self.field_area

        elif feature == VisualFeature.FIELD_AREA_POSITIONS:
            return self.field_area_positions

        elif feature == VisualFeature.COVERAGE:
            return self.converage

        else:
            raise ValueError("{} is a unknown visual feature".format(feature))

    def as_dict(self):
        """ordered dictionary with the most important feature"""
        rtn = [("Hash", self.oa.hash),
               ("Numerosity", self.numerosity),
               ("?", None),  # placeholder
               (VisualFeature.AV_PERIMETER.label(), self.average_perimeter),
               (VisualFeature.AV_SURFACE_AREA.label(), self.average_surface_area),
               (VisualFeature.TOTAL_PERIMETER.label(), self.total_perimeter),
               (VisualFeature.TOTAL_SURFACE_AREA.label(), self.total_surface_area),
               (VisualFeature.FIELD_AREA.label(), self.field_area),
               (VisualFeature.SPARSITY.label(), self.sparsity),
               (VisualFeature.COVERAGE.label(), self.converage),
               (VisualFeature.LOG_SIZE.label(), self.log_size),
               (VisualFeature.LOG_SPACING.label(), self.log_spacing)]

        if isinstance(self.oa, arrays.DotArray):
            rtn[2] = (VisualFeature.AV_DOT_DIAMETER.label(), self.average_dot_diameter)
        elif isinstance(self.oa, arrays.RectangleArray):
            rtn[2] = (VisualFeature.AV_RECT_SIZE.label(), self.average_rectangle_size)
        else:
            rtn.pop(2)
        return OrderedDict(rtn)

    def __str__(self):
        return self.as_text()

    def as_text(self, with_hash=True, extended_format=False, spacing_char="."):
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