# calculates visual features of a dot array/ dot cloud

from collections import OrderedDict
from enum import IntFlag, auto

import numpy as np
from scipy import spatial
from . import misc
from .geometry import cartesian2polar, polar2cartesian
from . import arrays


class VisualFeature(IntFlag):

    TOTAL_SURFACE_AREA = auto()
    ITEM_DIAMETER = auto()
    ITEM_SURFACE_AREA = auto()
    ITEM_PERIMETER = auto()
    TOTAL_PERIMETER = auto()
    RECT_SIZE = auto()
    SPARSITY = auto()
    FIELD_AREA = auto()
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
                        VisualFeature.ITEM_DIAMETER,
                        VisualFeature.ITEM_SURFACE_AREA,
                        VisualFeature.ITEM_PERIMETER,
                        VisualFeature.TOTAL_PERIMETER)

    def is_space_feature(self):
        return self in (VisualFeature.LOG_SPACING,
                        VisualFeature.SPARSITY,
                        VisualFeature.FIELD_AREA)

    def label(self):
        labels = {
            VisualFeature.LOG_SIZE: "Log Size",
            VisualFeature.TOTAL_SURFACE_AREA: "Total surface area",
            VisualFeature.ITEM_DIAMETER: "Mean item diameter",
            VisualFeature.ITEM_SURFACE_AREA: "Mean item surface area",
            VisualFeature.ITEM_PERIMETER: "Total perimeter",
            VisualFeature.TOTAL_PERIMETER: "Mean item perimeter",
            VisualFeature.RECT_SIZE: "Mean Rectangle Size",
            VisualFeature.LOG_SPACING: "Log Spacing",
            VisualFeature.SPARSITY: "Sparsity",
            VisualFeature.FIELD_AREA: "Field area",
            VisualFeature.COVERAGE: "Coverage"}
        return labels[self]


class ArrayFeatures(object):

    def __init__(self, object_array):
        # _lib or dot_cloud
        assert isinstance(object_array, (arrays.RectangleArray, arrays.DotArray))
        self.oa = object_array
        self._convex_hull = None

    def reset(self):
        """reset to enforce recalculation of certain parameter """
        self._convex_hull = None

    @property
    def convex_hull(self):
        if self._convex_hull is None:
            self._convex_hull = ConvexHull(self.oa)
        return self._convex_hull

    @property
    def mean_item_diameter(self):
        if not isinstance(self.oa, arrays.DotArray):
            return None
        return np.mean(self.oa.diameters)

    @property
    def mean_rectangle_size(self):
        if not isinstance(self.oa, arrays.RectangleArray):
            return None
        return np.mean(self.oa.sizes, axis=0)


    @property
    def total_surface_area(self):
        return np.sum(self.oa.surface_areas)

    @property
    def mean_item_surface_area(self):
        return np.mean(self.oa.surface_areas)

    @property
    def total_perimeter(self):
        return np.sum(self.oa.perimeter)

    @property
    def mean_item_perimeter(self):
        return np.mean(self.oa.perimeter)

    @property
    def field_area(self):
        return self.convex_hull.convex_hull_position.volume

    @property
    def numerosity(self):
        return len(self.oa._xy)

    @property
    def converage(self):
        """ percent coverage in the field area. It takes thus the item size
        into account. In contrast, the sparsity is only the ratio of field
        array and numerosity

        """

        try:
            return self.total_surface_area / self.field_area
        except:
            return None

    @property
    def log_size(self):
        return misc.log2(self.total_surface_area) + misc.log2(
            self.mean_item_surface_area)

    @property
    def log_spacing(self):
        return misc.log2(self.field_area) + misc.log2(self.sparsity)

    @property
    def sparsity(self):
        return self.field_area / self.numerosity

    @property
    def field_area_full(self):  # TODO not tested
        return self.convex_hull.outer_field_area

    def get(self, feature):
        """returns a feature"""

        assert isinstance(feature, VisualFeature)

       # Adapt
        if feature == VisualFeature.ITEM_DIAMETER:
            return self.mean_item_diameter

        elif feature == VisualFeature.ITEM_PERIMETER:
            return self.mean_item_perimeter

        elif feature == VisualFeature.TOTAL_PERIMETER:
            return self.total_perimeter

        elif feature == VisualFeature.ITEM_SURFACE_AREA:
            return self.mean_item_surface_area

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

        elif feature == VisualFeature.COVERAGE:
            return self.converage

        else:
            raise ValueError("{} is a unknown visual feature".format(feature))

    def as_dict(self):
        """ordered dictionary with the most important feature"""
        rtn = [("Hash", self.oa.hash),
               ("Numerosity", self.numerosity),
               (VisualFeature.TOTAL_SURFACE_AREA.label(), self.total_surface_area),
               (VisualFeature.ITEM_SURFACE_AREA.label(), self.mean_item_surface_area),
               ("?", None),  # placeholder
               (VisualFeature.ITEM_PERIMETER.label(), self.mean_item_perimeter),
               (VisualFeature.TOTAL_PERIMETER.label(), self.total_perimeter),
               (VisualFeature.FIELD_AREA.label(), self.field_area),
               (VisualFeature.SPARSITY.label(), self.sparsity),
               (VisualFeature.COVERAGE.label(), self.converage),
               (VisualFeature.LOG_SIZE.label(), self.log_size),
               (VisualFeature.LOG_SPACING.label(), self.log_spacing)]

        if isinstance(self.oa, arrays.DotArray):
            rtn[4] = (VisualFeature.ITEM_DIAMETER.label(), self.mean_item_diameter)
        elif isinstance(self.oa, arrays.RectangleArray):
            rtn[4] = (VisualFeature.RECT_SIZE.label(), self.mean_rectangle_size)
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

                    rtn += name + (spacing_char * (22 - len(name))) + (" " * (14 - len(value))) + value
        else:
            if with_hash:
                rtn = "ID: {} ".format(self.oa.hash)
            else:
                rtn = ""
            rtn += "N: {}, TSA: {}, ISA: {}, FA: {}, SPAR: {:.3f}, logSIZE: " \
                   "{:.2f}, logSPACE: {:.2f} COV: {:.2f}".format(
                self.numerosity,
                int(self.total_surface_area),
                int(self.mean_item_surface_area),
                int(self.field_area),
                self.sparsity,
                self.log_size,
                self.log_spacing,
                self.converage)
        return rtn


class ConvexHull(object):
    """convenient wrapper class for calculation of convex hulls
    """

    def __init__(self, object_array):
        assert isinstance(object_array, (arrays.RectangleArray, arrays.DotArray))

        self._xy = object_array.xy

        if isinstance(object_array, arrays.DotArray):
            # centered polar coordinates
            minmax = np.array((np.min(self._xy, axis=0), np.max(self._xy, axis=0)))
            center = np.reshape(minmax[1, :] - np.diff(minmax, axis=0) / 2, 2)  # center outer positions
            xy_centered = self._xy - center
            # outer positions
            polar_centered = cartesian2polar(xy_centered)
            polar_centered[:, 0] = polar_centered[:, 0] + (object_array.diameters / 2)
            self._outer_xy = polar2cartesian(polar_centered) + center

        elif isinstance(object_array, arrays.RectangleArray):
            # get all edges
            edges = []
            for rect in object_array.get():
                edges.extend([e.xy for e in rect.edges()])
            self._outer_xy = np.array(edges)

        self.convex_hull_position = spatial.ConvexHull(self._xy)
        self.convex_hull_outer = spatial.ConvexHull(self._outer_xy)

    @property
    def position_xy(self):
        return self._xy[self.convex_hull_position.vertices, :]

    @property
    def outer_xy(self):
        return self._outer_xy[self.convex_hull_outer.vertices, :]

    @property
    def outer_field_area(self):
        return self.convex_hull_outer.volume

    @property
    def position_field_area(self):
        return self.convex_hull_position.volume

